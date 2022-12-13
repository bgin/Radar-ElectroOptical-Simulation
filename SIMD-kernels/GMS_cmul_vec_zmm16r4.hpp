

#ifndef __GMS_CMUL_VEC_ZMM16R4_HPP__
#define __GMS_CMUL_VEC_ZMM16R4_HPP__

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

    const unsigned int GMS_CMUL_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CMUL_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CMUL_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CMUL_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CMUL_VEC_ZMM16R4_MAJOR+
      100U*GMS_CMUL_VEC_ZMM16R4_MINOR+
      10U*GMS_CMUL_VEC_ZMM16R4_MICRO;
    const char * const GMS_CMUL_VEC_ZMM16R4_CREATION_DATE = "12-12-2022 16:26 AM +00200 (MON 12 DEC 2022 GMT+2)";
    const char * const GMS_CMUL_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CMUL_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CMUL_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector multiplication operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"



namespace gms {


         namespace  math {


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cmul_zmm16r4_unroll_12x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) {return;}

                         register __m512 zmm0,zmm1,zmm2,zmm3;
                         register __m512 zmm4,zmm5,zmm6,zmm7;
                         register __m512 zmm8,zmm9,zmm10,zmm11;
                         register __m512 zmm12,zmm13,zmm14,zmm15;
                         register __m512 zmm16,zmm17,zmm18,zmm19;
                         register __m512 zmm20,zmm21,zmm22,zmm23;
                         register __m512 zmm24,zmm25,zmm26,zmm27;
                         register __m512 zmm28,zmm29;
                         int32_t i;

                         for(i = 0; (i+191) < n; i += 192) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm7 = _mm512_loadu_ps(&xre[i+96]);
                            zmm8 = _mm512_loadu_ps(&yre[i+96]);
                            zmm9 = _mm512_loadu_ps(&xim[i+96]);
                            zmm10= _mm512_loadu_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_storeu_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_loadu_ps(&xre[i+112]);
                            zmm14 = _mm512_loadu_ps(&yre[i+112]);
                            zmm15 = _mm512_loadu_ps(&xim[i+112]);
                            zmm16 = _mm512_loadu_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&zim[i+112], zmm18);
                            zmm19 = _mm512_loadu_ps(&xre[i+128]);
                            zmm20 = _mm512_loadu_ps(&yre[i+128]);
                            zmm21 = _mm512_loadu_ps(&xim[i+128]);
                            zmm22 = _mm512_loadu_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_storeu_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_add_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_storeu_ps(&zim[i+128], zmm24);
                            _mm_prefetch((const char *)&xre[i+196],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+196],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+196],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+196],_MM_HINT_T0);
                            zmm25 = _mm512_loadu_ps(&xre[i+144]);
                            zmm26 = _mm512_loadu_ps(&yre[i+144]);
                            zmm27 = _mm512_loadu_ps(&xim[i+144]);
                            zmm28 = _mm512_loadu_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_storeu_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_add_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_storeu_ps(&zim[i+144], zmm1);
                            zmm2 = _mm512_loadu_ps(&xre[i+160]);
                            zmm3 = _mm512_loadu_ps(&yre[i+160]);
                            zmm4 = _mm512_loadu_ps(&xim[i+160]);
                            zmm5 = _mm512_loadu_ps(&yim[i+160]);
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm2,zmm3),
                                                       _mm512_mul_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zre[i+160], zmm6);
                             zmm7  = _mm512_add_ps(_mm512_mul_ps(zmm4,zmm3),
                                                           _mm512_mul_ps(zmm2,zmm5));
                            _mm512_storeu_ps(&zim[i+160], zmm7);
                            zmm8 = _mm512_loadu_ps(&xre[i+176]);
                            zmm9 = _mm512_loadu_ps(&yre[i+176]);
                            zmm10= _mm512_loadu_ps(&xim[i+176]);
                            zmm11= _mm512_loadu_ps(&yim[i+176]);
                            zmm12= _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                       _mm512_mul_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+176], zmm12);
                             zmm13 = _mm512_add_ps(_mm512_mul_ps(zmm10,zmm8),
                                                           _mm512_mul_ps(zmm8,zmm11));
                            _mm512_storeu_ps(&zim[i+176], zmm13);

                         }

                      for(; (i+159) < n; i += 160) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_loadu_ps(&xre[i+96]);
                            zmm8 = _mm512_loadu_ps(&yre[i+96]);
                            zmm9 = _mm512_loadu_ps(&xim[i+96]);
                            zmm10= _mm512_loadu_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_storeu_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_loadu_ps(&xre[i+112]);
                            zmm14 = _mm512_loadu_ps(&yre[i+112]);
                            zmm15 = _mm512_loadu_ps(&xim[i+112]);
                            zmm16 = _mm512_loadu_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&zim[i+112], zmm18);
                            zmm19 = _mm512_loadu_ps(&xre[i+128]);
                            zmm20 = _mm512_loadu_ps(&yre[i+128]);
                            zmm21 = _mm512_loadu_ps(&xim[i+128]);
                            zmm22 = _mm512_loadu_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_storeu_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_add_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_storeu_ps(&zim[i+128], zmm24);
                            zmm25 = _mm512_loadu_ps(&xre[i+144]);
                            zmm26 = _mm512_loadu_ps(&yre[i+144]);
                            zmm27 = _mm512_loadu_ps(&xim[i+144]);
                            zmm28 = _mm512_loadu_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_storeu_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_add_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_storeu_ps(&zim[i+144], zmm1);

                      }

                    for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_loadu_ps(&xre[i+96]);
                            zmm8 = _mm512_loadu_ps(&yre[i+96]);
                            zmm9 = _mm512_loadu_ps(&xim[i+96]);
                            zmm10= _mm512_loadu_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_storeu_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_loadu_ps(&xre[i+112]);
                            zmm14 = _mm512_loadu_ps(&yre[i+112]);
                            zmm15 = _mm512_loadu_ps(&xim[i+112]);
                            zmm16 = _mm512_loadu_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_storeu_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_storeu_ps(&zim[i+112], zmm18);

                    }

                    for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_loadu_ps(&xre[i+80]);
                            zmm2 = _mm512_loadu_ps(&yre[i+80]);
                            zmm3 = _mm512_loadu_ps(&xim[i+80]);
                            zmm4 = _mm512_loadu_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_storeu_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                   }

                    for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);

                    }

                    for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                    }

                    for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);

                    }

                    for(; (i+0) < n; i += 1) {
                         zre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         zim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cmul_zmm16r4_unroll_12x_u(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yim,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) {return;}

                         register __m512 zmm0,zmm1,zmm2,zmm3;
                         register __m512 zmm4,zmm5,zmm6,zmm7;
                         register __m512 zmm8,zmm9,zmm10,zmm11;
                         register __m512 zmm12,zmm13,zmm14,zmm15;
                         register __m512 zmm16,zmm17,zmm18,zmm19;
                         register __m512 zmm20,zmm21,zmm22,zmm23;
                         register __m512 zmm24,zmm25,zmm26,zmm27;
                         register __m512 zmm28,zmm29;
                         int32_t i;

                         for(i = 0; (i+191) < n; i += 192) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&zim[i+32], zmm17);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_storeu_ps(&zim[i+80], zmm6);
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                            zmm7 = _mm512_load_ps(&xre[i+96]);
                            zmm8 = _mm512_load_ps(&yre[i+96]);
                            zmm9 = _mm512_load_ps(&xim[i+96]);
                            zmm10= _mm512_load_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_store_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+112]);
                            zmm14 = _mm512_load_ps(&yre[i+112]);
                            zmm15 = _mm512_load_ps(&xim[i+112]);
                            zmm16 = _mm512_load_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&zim[i+112], zmm18);
                            zmm19 = _mm512_load_ps(&xre[i+128]);
                            zmm20 = _mm512_load_ps(&yre[i+128]);
                            zmm21 = _mm512_load_ps(&xim[i+128]);
                            zmm22 = _mm512_load_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_store_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_add_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_store_ps(&zim[i+128], zmm24);
                            _mm_prefetch((const char *)&xre[i+196],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+196],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+196],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+196],_MM_HINT_T0);
                            zmm25 = _mm512_load_ps(&xre[i+144]);
                            zmm26 = _mm512_load_ps(&yre[i+144]);
                            zmm27 = _mm512_load_ps(&xim[i+144]);
                            zmm28 = _mm512_load_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_store_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_add_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_store_ps(&zim[i+144], zmm1);
                            zmm2 = _mm512_load_ps(&xre[i+160]);
                            zmm3 = _mm512_load_ps(&yre[i+160]);
                            zmm4 = _mm512_load_ps(&xim[i+160]);
                            zmm5 = _mm512_load_ps(&yim[i+160]);
                            zmm6 = _mm512_sub_ps(_mm512_mul_ps(zmm2,zmm3),
                                                       _mm512_mul_ps(zmm4,zmm5));
                            _mm512_store_ps(&zre[i+160], zmm6);
                             zmm7  = _mm512_add_ps(_mm512_mul_ps(zmm4,zmm3),
                                                           _mm512_mul_ps(zmm2,zmm5));
                            _mm512_store_ps(&zim[i+160], zmm7);
                            zmm8 = _mm512_load_ps(&xre[i+176]);
                            zmm9 = _mm512_load_ps(&yre[i+176]);
                            zmm10= _mm512_load_ps(&xim[i+176]);
                            zmm11= _mm512_load_ps(&yim[i+176]);
                            zmm12= _mm512_sub_ps(_mm512_mul_ps(zmm8,zmm9),
                                                       _mm512_mul_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+176], zmm12);
                             zmm13 = _mm512_add_ps(_mm512_mul_ps(zmm10,zmm8),
                                                           _mm512_mul_ps(zmm8,zmm11));
                            _mm512_store_ps(&zim[i+176], zmm13);

                         }

                      for(; (i+159) < n; i += 160) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_load_ps(&xre[i+96]);
                            zmm8 = _mm512_load_ps(&yre[i+96]);
                            zmm9 = _mm512_load_ps(&xim[i+96]);
                            zmm10= _mm512_load_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_store_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_store_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+112]);
                            zmm14 = _mm512_load_ps(&yre[i+112]);
                            zmm15 = _mm512_load_ps(&xim[i+112]);
                            zmm16 = _mm512_load_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&zim[i+112], zmm18);
                            zmm19 = _mm512_load_ps(&xre[i+128]);
                            zmm20 = _mm512_load_ps(&yre[i+128]);
                            zmm21 = _mm512_load_ps(&xim[i+128]);
                            zmm22 = _mm512_load_ps(&yim[i+128]);
                            zmm23 = _mm512_sub_ps(_mm512_mul_ps(zmm19,zmm20),
                                                       _mm512_mul_ps(zmm21,zmm22));
                            _mm512_store_ps(&zre[i+128], zmm23);
                             zmm24  = _mm512_add_ps(_mm512_mul_ps(zmm21,zmm20),
                                                           _mm512_mul_ps(zmm19,zmm22));
                            _mm512_store_ps(&zim[i+128], zmm24);
                            zmm25 = _mm512_load_ps(&xre[i+144]);
                            zmm26 = _mm512_load_ps(&yre[i+144]);
                            zmm27 = _mm512_load_ps(&xim[i+144]);
                            zmm28 = _mm512_load_ps(&yim[i+144]);
                            zmm29 = _mm512_sub_ps(_mm512_mul_ps(zmm25,zmm26),
                                                       _mm512_mul_ps(zmm27,zmm28));
                            _mm512_store_ps(&zre[i+144], zmm29);
                             zmm1  = _mm512_add_ps(_mm512_mul_ps(zmm27,zmm26),
                                                           _mm512_mul_ps(zmm25,zmm28));
                            _mm512_store_ps(&zim[i+144], zmm1);

                      }

                    for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_store_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&zim[i+80], zmm6);
                            zmm7 = _mm512_load_ps(&xre[i+96]);
                            zmm8 = _mm512_load_ps(&yre[i+96]);
                            zmm9 = _mm512_load_ps(&xim[i+96]);
                            zmm10= _mm512_load_ps(&yim[i+96]);
                            zmm11= _mm512_sub_ps(_mm512_mul_ps(zmm7,zmm8),
                                                       _mm512_mul_ps(zmm9,zmm10));
                            _mm512_store_ps(&zre[i+96], zmm11);
                             zmm12  = _mm512_add_ps(_mm512_mul_ps(zmm9,zmm8),
                                                           _mm512_mul_ps(zmm7,zmm10));
                            _mm512_storeu_ps(&zim[i+96], zmm12);
                            zmm13 = _mm512_load_ps(&xre[i+112]);
                            zmm14 = _mm512_load_ps(&yre[i+112]);
                            zmm15 = _mm512_load_ps(&xim[i+112]);
                            zmm16 = _mm512_load_ps(&yim[i+112]);
                            zmm17 = _mm512_sub_ps(_mm512_mul_ps(zmm13,zmm14),
                                                       _mm512_mul_ps(zmm15,zmm16));
                            _mm512_store_ps(&zre[i+112], zmm17;
                             zmm18  = _mm512_add_ps(_mm512_mul_ps(zmm15,zmm14),
                                                           _mm512_mul_ps(zmm13,zmm16));
                            _mm512_store_ps(&zim[i+112], zmm18);

                    }

                    for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_store_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_storeu_ps(&zim[i+48], zmm23);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_sub_ps(_mm512_mul_ps(zmm24,zmm25),
                                                                        _mm512_mul_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+64], zmm28);
                             zmm29  = _mm512_add_ps(_mm512_mul_ps(zmm26,zmm25),
                                                                        _mm512_mul_ps(zmm24,zmm27));
                            _mm512_storeu_ps(&zim[i+64], zmm29);
                            zmm1 = _mm512_load_ps(&xre[i+80]);
                            zmm2 = _mm512_load_ps(&yre[i+80]);
                            zmm3 = _mm512_load_ps(&xim[i+80]);
                            zmm4 = _mm512_load_ps(&yim[i+80]);
                            zmm5 = _mm512_sub_ps(_mm512_mul_ps(zmm1,zmm2),
                                                       _mm512_mul_ps(zmm3,zmm4));
                            _mm512_store_ps(&zre[i+80], zmm5);
                             zmm6  = _mm512_add_ps(_mm512_mul_ps(zmm3,zmm2),
                                                           _mm512_mul_ps(zmm1,zmm4));
                            _mm512_store_ps(&zim[i+80], zmm6);
                   }

                    for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_storeu_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_storeu_ps(&zim[i+16], zmm11);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16 = _mm512_sub_ps(_mm512_mul_ps(zmm12,zmm13),
                                                                        _mm512_mul_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+32], zmm16);
                             zmm17  = _mm512_add_ps(_mm512_mul_ps(zmm14,zmm13),
                                                                        _mm512_mul_ps(zmm12,zmm15));
                            _mm512_storeu_ps(&zim[i+32], zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_sub_ps(_mm512_mul_ps(zmm18,zmm19),
                                                                        _mm512_mul_ps(zmm20,zmm21));
                            _mm512_store_ps(&zre[i+48], zmm22);
                             zmm23  = _mm512_add_ps(_mm512_mul_ps(zmm20,zmm19),
                                                                        _mm512_mul_ps(zmm18,zmm21));
                            _mm512_store_ps(&zim[i+48], zmm23);

                    }

                    for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             zmm10 = _mm512_sub_ps(_mm512_mul_ps(zmm6,zmm7),
                                                                        _mm512_mul_ps(zmm8,zmm9));
                            _mm512_store_ps(&zre[i+16], zmm10);
                            zmm11 = _mm512_add_ps(_mm512_mul_ps(zmm8,zmm7),
                                                                        _mm512_mul_ps(zmm6,zmm9));
                            _mm512_store_ps(&zim[i+16], zmm11);
                    }

                    for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+0], zmm4);
                             zmm5  = _mm512_add_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                            _mm512_store_ps(&zim[i+0], zmm5);

                    }

                    for(; (i+0) < n; i += 1) {
                         zre[i] = (xre[i] * yre[i]) - (xim[i] * yim[i]);
		         zim[i] = (xim[i] * yre[i]) + (xre[i] * yim[i]);
                    }
               }



      
        } // math


} // gms















#endif /*__GMS_CMUL_VEC_ZMM16R4_HPP__*/
