


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
#include "GMS_cabs_vec_zmm16r4.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"




	           
	                 void gms::math::cabsv_zmm16r4_unroll_16x_u(const float * __restrict re,
                                                   const float * __restrict im,
                                                   float * __restrict  cabs,
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
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_loadu_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_loadu_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_storeu_ps(&cabs[i+80],zmm29);
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm30  = _mm512_loadu_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_loadu_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_storeu_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_loadu_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_loadu_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_storeu_ps(&cabs[i+112],zmm7);
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&re[i+128]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&im[i+128]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = xsqrtf(_mm512_add_ps(zmm9,zmm11));
                              _mm512_storeu_ps(&cabs[i+128],zmm12);
                              zmm13 = _mm512_loadu_ps(&re[i+144]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = _mm512_loadu_ps(&im[i+144]);
                              zmm16 = _mm512_mul_ps(zmm15,zmm15);
                              zmm17 = xsqrtf(_mm512_add_ps(zmm14,zmm16));
                              _mm512_storeu_ps(&cabs[i+144],zmm17);
                              _mm_prefetch((const char*)&re[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+192],_MM_HINT_T0);
                              zmm18 = _mm512_loadu_ps(&re[i+160]);
                              zmm19 = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_loadu_ps(&im[i+160]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = xsqrtf(_mm512_add_ps(zmm19,zmm21));
                              _mm512_storeu_ps(&cabs[i+160],zmm22);
                              zmm23 = _mm512_loadu_ps(&re[i+176]);
                              zmm24 = _mm512_mul_ps(zmm23,zmm23);
                              zmm25 = _mm512_loadu_ps(&im[i+176]);
                              zmm26 = _mm512_mul_ps(zmm25,zmm25);
                              zmm27 = xsqrtf(_mm512_add_ps(zmm24,zmm26));
                              _mm512_storeu_ps(&cabs[i+176],zmm27);
                              _mm_prefetch((const char*)&re[i+224],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+224],_MM_HINT_T0);
                              zmm28 = _mm512_loadu_ps(&re[i+192]);
                              zmm29 = _mm512_mul_ps(zmm28,zmm28);
                              zmm30 = _mm512_loadu_ps(&im[i+192]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = xsqrtf(_mm512_add_ps(zmm29,zmm31));
                              _mm512_storeu_ps(&cabs[i+192],zmm0);
                              zmm1  = _mm512_loadu_ps(&re[i+208]);
                              zmm2  = _mm512_mul_ps(zmm1,zmm1);
                              zmm3  = _mm512_loadu_ps(&im[i+208]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = xsqrtf(_mm512_add_ps(zmm2,zmm4));
                              _mm512_storeu_ps(&cabs[i+208],zmm5);
                              _mm_prefetch((const char*)&re[i+256],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+256],_MM_HINT_T0);
                              zmm6  = _mm512_loadu_ps(&re[i+224]);
                              zmm7  = _mm512_mul_ps(zmm6,zmm6);
                              zmm8  = _mm512_loadu_ps(&im[i+224]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = xsqrtf(_mm512_add_ps(zmm7,zmm9));
                              _mm512_storeu_ps(&cabs[i+224],zmm10);
                              zmm11 = _mm512_loadu_ps(&re[i+240]);
                              zmm12 = _mm512_mul_ps(zmm11,zmm11);
                              zmm13 = _mm512_loadu_ps(&im[i+240]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = xsqrtf(_mm512_add_ps(zmm12,zmm14));
                              _mm512_storeu_ps(&cabs[i+240],zmm15);
                        }

                         for(; (i+191) < n; i += 192) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_loadu_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_loadu_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_storeu_ps(&cabs[i+80],zmm29);
                              zmm30  = _mm512_loadu_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_loadu_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_storeu_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_loadu_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_loadu_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_storeu_ps(&cabs[i+112],zmm7);
                              zmm8  = _mm512_loadu_ps(&re[i+128]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&im[i+128]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = xsqrtf(_mm512_add_ps(zmm9,zmm11));
                              _mm512_storeu_ps(&cabs[i+128],zmm12);
                              zmm13 = _mm512_loadu_ps(&re[i+144]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = _mm512_loadu_ps(&im[i+144]);
                              zmm16 = _mm512_mul_ps(zmm15,zmm15);
                              zmm17 = xsqrtf(_mm512_add_ps(zmm14,zmm16));
                              _mm512_storeu_ps(&cabs[i+144],zmm17);
                              zmm18 = _mm512_loadu_ps(&re[i+160]);
                              zmm19 = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_loadu_ps(&im[i+160]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = xsqrtf(_mm512_add_ps(zmm19,zmm21));
                              _mm512_storeu_ps(&cabs[i+160],zmm22);
                              zmm23 = _mm512_loadu_ps(&re[i+176]);
                              zmm24 = _mm512_mul_ps(zmm23,zmm23);
                              zmm25 = _mm512_loadu_ps(&im[i+176]);
                              zmm26 = _mm512_mul_ps(zmm25,zmm25);
                              zmm27 = xsqrtf(_mm512_add_ps(zmm24,zmm26));
                              _mm512_storeu_ps(&cabs[i+176],zmm27);
                        }

                         for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_loadu_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_loadu_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_storeu_ps(&cabs[i+80],zmm29);
                              zmm30  = _mm512_loadu_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_loadu_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_storeu_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_loadu_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_loadu_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_storeu_ps(&cabs[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                        }

                         for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float re2 = xre*xre;
                              const float xim = im[i];
                              const float im2 = xim*xim;
                              cabs[i]         = ceph_sqrtf(re2+im2);
                        }


                
                   void gms::math::cabsv_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict  __ATTR_ALIGN__(64) cabs,
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
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_store_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_load_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_load_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_store_ps(&cabs[i+80],zmm29);
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm30  = _mm512_load_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_load_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_store_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_load_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_load_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_store_ps(&cabs[i+112],zmm7);
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&re[i+128]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&im[i+128]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = xsqrtf(_mm512_add_ps(zmm9,zmm11));
                              _mm512_store_ps(&cabs[i+128],zmm12);
                              zmm13 = _mm512_load_ps(&re[i+144]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = _mm512_load_ps(&im[i+144]);
                              zmm16 = _mm512_mul_ps(zmm15,zmm15);
                              zmm17 = xsqrtf(_mm512_add_ps(zmm14,zmm16));
                              _mm512_store_ps(&cabs[i+144],zmm17);
                              _mm_prefetch((const char*)&re[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+192],_MM_HINT_T0);
                              zmm18 = _mm512_load_ps(&re[i+160]);
                              zmm19 = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_load_ps(&im[i+160]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = xsqrtf(_mm512_add_ps(zmm19,zmm21));
                              _mm512_store_ps(&cabs[i+160],zmm22);
                              zmm23 = _mm512_load_ps(&re[i+176]);
                              zmm24 = _mm512_mul_ps(zmm23,zmm23);
                              zmm25 = _mm512_load_ps(&im[i+176]);
                              zmm26 = _mm512_mul_ps(zmm25,zmm25);
                              zmm27 = xsqrtf(_mm512_add_ps(zmm24,zmm26));
                              _mm512_store_ps(&cabs[i+176],zmm27);
                              _mm_prefetch((const char*)&re[i+224],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+224],_MM_HINT_T0);
                              zmm28 = _mm512_load_ps(&re[i+192]);
                              zmm29 = _mm512_mul_ps(zmm28,zmm28);
                              zmm30 = _mm512_load_ps(&im[i+192]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = xsqrtf(_mm512_add_ps(zmm29,zmm31));
                              _mm512_store_ps(&cabs[i+192],zmm0);
                              zmm1  = _mm512_load_ps(&re[i+208]);
                              zmm2  = _mm512_mul_ps(zmm1,zmm1);
                              zmm3  = _mm512_load_ps(&im[i+208]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = xsqrtf(_mm512_add_ps(zmm2,zmm4));
                              _mm512_store_ps(&cabs[i+208],zmm5);
                              _mm_prefetch((const char*)&re[i+256],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+256],_MM_HINT_T0);
                              zmm6  = _mm512_load_ps(&re[i+224]);
                              zmm7  = _mm512_mul_ps(zmm6,zmm6);
                              zmm8  = _mm512_load_ps(&im[i+224]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = xsqrtf(_mm512_add_ps(zmm7,zmm9));
                              _mm512_store_ps(&cabs[i+224],zmm10);
                              zmm11 = _mm512_load_ps(&re[i+240]);
                              zmm12 = _mm512_mul_ps(zmm11,zmm11);
                              zmm13 = _mm512_load_ps(&im[i+240]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = xsqrtf(_mm512_add_ps(zmm12,zmm14));
                              _mm512_store_ps(&cabs[i+240],zmm15);
                        }

                         for(; (i+191) < n; i += 192) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_store_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_load_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_load_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_store_ps(&cabs[i+80],zmm29);
                              zmm30  = _mm512_load_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_load_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_store_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_load_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_load_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_store_ps(&cabs[i+112],zmm7);
                              zmm8  = _mm512_load_ps(&re[i+128]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&im[i+128]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = xsqrtf(_mm512_add_ps(zmm9,zmm11));
                              _mm512_store_ps(&cabs[i+128],zmm12);
                              zmm13 = _mm512_load_ps(&re[i+144]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = _mm512_load_ps(&im[i+144]);
                              zmm16 = _mm512_mul_ps(zmm15,zmm15);
                              zmm17 = xsqrtf(_mm512_add_ps(zmm14,zmm16));
                              _mm512_store_ps(&cabs[i+144],zmm17);
                              zmm18 = _mm512_load_ps(&re[i+160]);
                              zmm19 = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_load_ps(&im[i+160]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = xsqrtf(_mm512_add_ps(zmm19,zmm21));
                              _mm512_store_ps(&cabs[i+160],zmm22);
                              zmm23 = _mm512_load_ps(&re[i+176]);
                              zmm24 = _mm512_mul_ps(zmm23,zmm23);
                              zmm25 = _mm512_load_ps(&im[i+176]);
                              zmm26 = _mm512_mul_ps(zmm25,zmm25);
                              zmm27 = xsqrtf(_mm512_add_ps(zmm24,zmm26));
                              _mm512_store_ps(&cabs[i+176],zmm27);
                        }

                         for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_load_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_load_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_store_ps(&cabs[i+80],zmm29);
                              zmm30  = _mm512_load_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_load_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_store_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_load_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_load_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_store_ps(&cabs[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                        }

                         for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float re2 = xre*xre;
                              const float xim = im[i];
                              const float im2 = xim*xim;
                              cabs[i]         = ceph_sqrtf(re2+im2);
                        }

              }


                  
                   void gms::math::cabsv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict  __ATTR_ALIGN__(64) cabs,
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
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_store_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_load_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_load_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_store_ps(&cabs[i+80],zmm29);
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm30  = _mm512_load_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_load_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_store_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_load_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_load_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_store_ps(&cabs[i+112],zmm7);
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&re[i+128]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&im[i+128]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = xsqrtf(_mm512_add_ps(zmm9,zmm11));
                              _mm512_store_ps(&cabs[i+128],zmm12);
                              zmm13 = _mm512_load_ps(&re[i+144]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = _mm512_load_ps(&im[i+144]);
                              zmm16 = _mm512_mul_ps(zmm15,zmm15);
                              zmm17 = xsqrtf(_mm512_add_ps(zmm14,zmm16));
                              _mm512_store_ps(&cabs[i+144],zmm17);
                             
                        }

                       
                         for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_load_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_load_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_store_ps(&cabs[i+80],zmm29);
                              zmm30  = _mm512_load_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_load_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_store_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_load_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_load_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_store_ps(&cabs[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                        }

                         for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float re2 = xre*xre;
                              const float xim = im[i];
                              const float im2 = xim*xim;
                              cabs[i]         = ceph_sqrtf(re2+im2);
                        }

              }


                 
                   void gms::math::cabsv_zmm16r4_unroll_10x_u(const float * __restrict re,
                                                   const float * __restrict im,
                                                   float * __restrict  cabs,
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
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_loadu_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_loadu_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_storeu_ps(&cabs[i+80],zmm29);
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm30  = _mm512_loadu_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_loadu_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_storeu_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_loadu_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_loadu_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_storeu_ps(&cabs[i+112],zmm7);
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&re[i+128]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&im[i+128]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = xsqrtf(_mm512_add_ps(zmm9,zmm11));
                              _mm512_storeu_ps(&cabs[i+128],zmm12);
                              zmm13 = _mm512_loadu_ps(&re[i+144]);
                              zmm14 = _mm512_mul_ps(zmm13,zmm13);
                              zmm15 = _mm512_loadu_ps(&im[i+144]);
                              zmm16 = _mm512_mul_ps(zmm15,zmm15);
                              zmm17 = xsqrtf(_mm512_add_ps(zmm14,zmm16));
                              _mm512_storeu_ps(&cabs[i+144],zmm17);
                           
                        }

                       
                         for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_loadu_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_loadu_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_storeu_ps(&cabs[i+80],zmm29);
                              zmm30  = _mm512_loadu_ps(&re[i+96]);
                              zmm31  = _mm512_mul_ps(zmm30,zmm30);
                              zmm0   = _mm512_loadu_ps(&im[i+96]);
                              zmm1   = _mm512_mul_ps(zmm0,zmm0);
                              zmm2   = xsqrtf(_mm512_add_ps(zmm31,zmm1));
                              _mm512_storeu_ps(&cabs[i+96],zmm2);
                              zmm3  = _mm512_loadu_ps(&re[i+112]);
                              zmm4  = _mm512_mul_ps(zmm3,zmm3);
                              zmm5  = _mm512_loadu_ps(&im[i+112]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = xsqrtf(_mm512_add_ps(zmm4,zmm6));
                              _mm512_storeu_ps(&cabs[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                        }

                         for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float re2 = xre*xre;
                              const float xim = im[i];
                              const float im2 = xim*xim;
                              cabs[i]         = ceph_sqrtf(re2+im2);
                        }

                }


                 
                   void gms::math::cabsv_zmm16r4_unroll_6x_u( const float * __restrict re,
                                                   const float * __restrict im,
                                                   float * __restrict  cabs,
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
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_loadu_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_loadu_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_storeu_ps(&cabs[i+80],zmm29);
                             
                        }

                       
                        

                         for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_loadu_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_loadu_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_storeu_ps(&cabs[i+64],zmm23);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_loadu_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_loadu_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_loadu_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_loadu_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_storeu_ps(&cabs[i+48],zmm19);
                        }

                         for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_loadu_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_loadu_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float re2 = xre*xre;
                              const float xim = im[i];
                              const float im2 = xim*xim;
                              cabs[i]         = ceph_sqrtf(re2+im2);
                        }

                }


                
                   void gms::math::cabsv_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict  __ATTR_ALIGN__(64) cabs,
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
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_store_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                              zmm25  = _mm512_load_ps(&re[i+80]);
                              zmm26  = _mm512_mul_ps(zmm25,zmm25);
                              zmm27  = _mm512_load_ps(&im[i+80]);
                              zmm28  = _mm512_mul_ps(zmm27,zmm27);
                              zmm29  = xsqrtf(_mm512_add_ps(zmm26,zmm28));
                              _mm512_store_ps(&cabs[i+80],zmm29);
                                                          
                        }

                       
                      

                         for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                              zmm20  = _mm512_load_ps(&re[i+64]);
                              zmm21  = _mm512_mul_ps(zmm21,zmm21);
                              zmm22  = _mm512_load_ps(&im[i+64]);
                              zmm23  = _mm512_mul_ps(zmm22,zmm22);
                              zmm24  = xsqrtf(_mm512_add_ps(zmm21,zmm23));
                              _mm512_store_ps(&cabs[i+64],zmm23);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_storeu_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_storeu_ps(&cabs[i+16],zmm9);
                              zmm10  = _mm512_load_ps(&re[i+32]);
                              zmm11  = _mm512_mul_ps(zmm10,zmm10);
                              zmm12  = _mm512_load_ps(&im[i+32]);
                              zmm13  = _mm512_mul_ps(zmm12,zmm12);
                              zmm14  = xsqrtf(_mm512_add_ps(zmm11,zmm13));
                              _mm512_storeu_ps(&cabs[i+32],zmm14);
                              zmm15  = _mm512_load_ps(&re[i+48]);
                              zmm16  = _mm512_mul_ps(zmm15,zmm15);
                              zmm17  = _mm512_load_ps(&im[i+48]);
                              zmm18  = _mm512_mul_ps(zmm17,zmm17);
                              zmm19  = xsqrtf(_mm512_add_ps(zmm16,zmm18));
                              _mm512_store_ps(&cabs[i+48],zmm19);
                        }

                         for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                              zmm5  = _mm512_load_ps(&re[i+16]);
                              zmm6  = _mm512_mul_ps(zmm5,zmm5);
                              zmm7  = _mm512_load_ps(&im[i+16]);
                              zmm8  = _mm512_mul_ps(zmm7,zmm7);
                              zmm9  = xsqrtf(_mm512_add_ps(zmm6,zmm8));
                              _mm512_store_ps(&cabs[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&im[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                              _mm512_store_ps(&cabs[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float re2 = xre*xre;
                              const float xim = im[i];
                              const float im2 = xim*xim;
                              cabs[i]         = ceph_sqrtf(re2+im2);
                        }

              }



                 



