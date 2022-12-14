

#ifndef __GMS_CADD_VEC_ZMM16R4_HPP__
#define __GMS_CADD_VEC_ZMM16R4_HPP__ 271120220939


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

    const unsigned int GMS_CADD_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CADD_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CADD_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CADD_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CADD_VEC_ZMM16R4_MAJOR+
      100U*GMS_CADD_VEC_ZMM16R4_MINOR+
      10U*GMS_CADD_VEC_ZMM16R4_MICRO;
    const char * const GMS_CADD_VEC_ZMM16R4_CREATION_DATE = "27-11-2022 09:39 AM +00200 (SUN 27 NOV 2022 GMT+2)";
    const char * const GMS_CADD_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CADD_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CADD_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector addition operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

         namespace math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_16x_u(const float * __restrict  xre,
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
                        register __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));
                            zmm0 = _mm512_loadu_ps(&xre[i+128]); 
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_add_ps(zmm0,zmm1));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_add_ps(zmm2,zmm3));
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            zmm4 = _mm512_loadu_ps(&xre[i+144]); 
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_add_ps(zmm4,zmm5));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_add_ps(zmm6,zmm7));
                            zmm8 = _mm512_loadu_ps(&xre[i+160]); 
                            zmm9 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+160]);
                            zmm11 = _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+176]); 
                            zmm13 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+176]);
                            zmm15 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+192]); 
                            zmm17 = _mm512_loadu_ps(&yre[i+192]);
                            _mm512_storeu_ps(&zre[i+192], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+192]);
                            zmm19 = _mm512_loadu_ps(&yim[i+192]);
                            _mm512_storeu_ps(&zim[i+192], _mm512_add_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+255],_MM_HINT_T0);
                            zmm20 = _mm512_loadu_ps(&xre[i+208]); 
                            zmm21 = _mm512_loadu_ps(&yre[i+208]);
                            _mm512_storeu_ps(&zre[i+208], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+208]);
                            zmm23 = _mm512_loadu_ps(&yim[i+208]);
                            _mm512_storeu_ps(&zim[i+208], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+224]); 
                            zmm25 = _mm512_loadu_ps(&yre[i+224]);
                            _mm512_storeu_ps(&zre[i+224], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+224]);
                            zmm27 = _mm512_loadu_ps(&yim[i+224]);
                            _mm512_storeu_ps(&zim[i+224], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+240]); 
                            zmm29 = _mm512_loadu_ps(&yre[i+240]);
                            _mm512_storeu_ps(&zre[i+240], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+240]);
                            zmm31 = _mm512_loadu_ps(&yim[i+240]);
                            _mm512_storeu_ps(&zim[i+240], _mm512_add_ps(zmm30,zmm31));
                        
                        }

                       for(; (i+191) < n; i += 192) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));
                            zmm0 = _mm512_loadu_ps(&xre[i+128]); 
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_add_ps(zmm0,zmm1));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_add_ps(zmm2,zmm3));
                            zmm4 = _mm512_loadu_ps(&xre[i+144]); 
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_add_ps(zmm4,zmm5));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_add_ps(zmm6,zmm7));
                            zmm8 = _mm512_loadu_ps(&xre[i+160]); 
                            zmm9 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+160]);
                            zmm11 = _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+176]); 
                            zmm13 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+176]);
                            zmm15 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_add_ps(zmm14,zmm15));

                       }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));

                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));

                       }

                      for(; (i+15) < n; i += 16) {
                           zmm0  = _mm512_loadu_ps(&xre[i+0]);
                           zmm1  = _mm512_loadu_ps(&yre[i+0]);
                           _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                           zmm2  = _mm512_loadu_ps(&xim[i+0]);
                           zmm3  = _mm512_loadu_ps(&yim[i+0]);
                           _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] + yre[i];
                            zim[i] = xim[i] + yim[i];
                      }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_16x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
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
                        register __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));
                            zmm0 = _mm512_load_ps(&xre[i+128]); 
                            zmm1 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_add_ps(zmm0,zmm1));
                            zmm2 = _mm512_load_ps(&xim[i+128]);
                            zmm3 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_add_ps(zmm2,zmm3));
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            zmm4 = _mm512_load_ps(&xre[i+144]); 
                            zmm5 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_add_ps(zmm4,zmm5));
                            zmm6 = _mm512_load_ps(&xim[i+144]);
                            zmm7 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_add_ps(zmm6,zmm7));
                            zmm8 = _mm512_load_ps(&xre[i+160]); 
                            zmm9 = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+160]);
                            zmm11 = _mm512_load_ps(&yim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+176]); 
                            zmm13 = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+176]);
                            zmm15 = _mm512_load_ps(&yim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+192]); 
                            zmm17 = _mm512_load_ps(&yre[i+192]);
                            _mm512_store_ps(&zre[i+192], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+192]);
                            zmm19 = _mm512_load_ps(&yim[i+192]);
                            _mm512_store_ps(&zim[i+192], _mm512_add_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+255],_MM_HINT_T0);
                            zmm20 = _mm512_load_ps(&xre[i+208]); 
                            zmm21 = _mm512_load_ps(&yre[i+208]);
                            _mm512_store_ps(&zre[i+208], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+208]);
                            zmm23 = _mm512_load_ps(&yim[i+208]);
                            _mm512_store_ps(&zim[i+208], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+224]); 
                            zmm25 = _mm512_load_ps(&yre[i+224]);
                            _mm512_store_ps(&zre[i+224], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+224]);
                            zmm27 = _mm512_load_ps(&yim[i+224]);
                            _mm512_store_ps(&zim[i+224], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+240]); 
                            zmm29 = _mm512_load_ps(&yre[i+240]);
                            _mm512_store_ps(&zre[i+240], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+240]);
                            zmm31 = _mm512_load_ps(&yim[i+240]);
                            _mm512_store_ps(&zim[i+240], _mm512_add_ps(zmm30,zmm31));
                        
                        }

                       for(; (i+191) < n; i += 192) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));
                            zmm0 = _mm512_load_ps(&xre[i+128]); 
                            zmm1 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_add_ps(zmm0,zmm1));
                            zmm2 = _mm512_load_ps(&xim[i+128]);
                            zmm3 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_add_ps(zmm2,zmm3));
                            zmm4 = _mm512_load_ps(&xre[i+144]); 
                            zmm5 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_add_ps(zmm4,zmm5));
                            zmm6 = _mm512_load_ps(&xim[i+144]);
                            zmm7 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_add_ps(zmm6,zmm7));
                            zmm8 = _mm512_load_ps(&xre[i+160]); 
                            zmm9 = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+160]);
                            zmm11 = _mm512_load_ps(&yim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+176]); 
                            zmm13 = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+176]);
                            zmm15 = _mm512_load_ps(&yim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_add_ps(zmm14,zmm15));
                       }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));

                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));

                       }

                      for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));

                      }

                      for(; (i+0) < n; i += 1) {
                            zre[i] = xre[i] + yre[i];
                            zim[i] = xim[i] + yim[i];
                      }
                  }

                    

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_10x_u(const float * __restrict  xre,
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
                        register __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;

                        for(i = 0; (i+159) < n; i += 160) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));
                            zmm0 = _mm512_loadu_ps(&xre[i+128]); 
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_add_ps(zmm0,zmm1));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_add_ps(zmm2,zmm3));
                            zmm4 = _mm512_loadu_ps(&xre[i+144]); 
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_add_ps(zmm4,zmm5));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_add_ps(zmm6,zmm7));
                        }


                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));

                       }

                       for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                                                 
                       }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                      }

                      for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                      }

                     for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                     }

                     for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] + yre[i];
                            zim[i] = xim[i] + yim[i];
                      }

                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_10x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
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
                        register __m512 zmm28,zmm29,zmm30,zmm31;
                        int32_t i;

                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                             zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                             zmm22 = _mm512_load_ps(&xim[i+80]);
                             zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                             zmm24 = _mm512_load_ps(&xre[i+96]);
                             zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                             zmm26 = _mm512_load_ps(&xim[i+96]);
                             zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                             zmm28 = _mm512_load_ps(&xre[i+112]);
                             zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                             zmm30 = _mm512_load_ps(&xim[i+112]);
                             zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));
                             zmm32 = _mm512_load_ps(&xre[i+128]);
                             zmm33 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_add_ps(zmm32,zmm33));
                             zmm34 = _mm512_load_ps(&xim[i+128]);
                             zmm35 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_add_ps(zmm34,zmm35));
                             zmm36 = _mm512_load_ps(&xre[i+144]); 
                             zmm37 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_add_ps(zmm36,zmm37));
                             zmm38 = _mm512_load_ps(&xim[i+144]);
                             zmm39 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_add_ps(zmm38,zmm39));

                        }

                       for(; (i+127) < n; i += 128) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                            zmm24 = _mm512_load_ps(&xre[i+96]);
                            zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm24,zmm25));
                            zmm26 = _mm512_load_ps(&xim[i+96]);
                            zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm26,zmm27));
                            zmm28 = _mm512_load_ps(&xre[i+112]);
                            zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm28,zmm29));
                            zmm30 = _mm512_load_ps(&xim[i+112]);
                            zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm30,zmm31));

                       }

                       for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                             zmm4  = _mm512_load_ps(&xre[i+16]);
                             zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                             zmm8  = _mm512_load_ps(&xre[i+32]);
                             zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                             zmm10 = _mm512_load_ps(&xim[i+32]);
                             zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                             zmm12 = _mm512_load_ps(&xre[i+48]);
                             zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                             zmm14 = _mm512_load_ps(&xim[i+48]);
                             zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                             zmm16 = _mm512_load_ps(&xre[i+64]);
                             zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                             zmm18 = _mm512_load_ps(&xim[i+64]);
                             zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));

                       }

                       for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                             zmm4  = _mm512_load_ps(&xre[i+16]);
                             zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                             zmm8  = _mm512_load_ps(&xre[i+32]);
                             zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                             zmm10 = _mm512_load_ps(&xim[i+32]);
                             zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                             zmm12 = _mm512_load_ps(&xre[i+48]);
                             zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                             zmm14 = _mm512_load_ps(&xim[i+48]);
                             zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));

                      }

                      for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                             zmm4  = _mm512_load_ps(&xre[i+16]);
                             zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));

                      }

                     for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));

                     }

                     for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] + yre[i];
                            zim[i] = xim[i] + yim[i];
                      }

                  }
 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_6x_u(const float * __restrict  xre,
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
                        int32_t i;

                        for(i = 0; (i+95) < n; i += 96) {
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                             zmm4  = _mm512_loadu_ps(&xre[i+16]);
                             zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                             zmm8  = _mm512_loadu_ps(&xre[i+32]);
                             zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                             zmm10 = _mm512_loadu_ps(&xim[i+32]);
                             zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                             zmm12 = _mm512_loadu_ps(&xre[i+48]);
                             zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                             zmm14 = _mm512_loadu_ps(&xim[i+48]);
                             zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0); 
                             zmm16 = _mm512_loadu_ps(&xre[i+64]);
                             zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                             zmm18 = _mm512_loadu_ps(&xim[i+64]);
                             zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                             zmm20 = _mm512_loadu_ps(&xre[i+80]);
                             zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                             zmm22 = _mm512_loadu_ps(&xim[i+80]);
                             zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));
                        }

                       for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                             zmm4  = _mm512_loadu_ps(&xre[i+16]);
                             zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                             zmm8  = _mm512_loadu_ps(&xre[i+32]);
                             zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                             zmm10 = _mm512_loadu_ps(&xim[i+32]);
                             zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                             zmm12 = _mm512_loadu_ps(&xre[i+48]);
                             zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                             zmm14 = _mm512_loadu_ps(&xim[i+48]);
                             zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                       }

                       for(; (i+31) < n; i += 32) {

                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                             zmm4  = _mm512_loadu_ps(&xre[i+16]);
                             zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                       }

                       for(; (i+15) < n; i += 16) {

                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));

                       }

                        for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] + yre[i];
                            zim[i] = xim[i] + yim[i];
                      }

                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_6x_a( const float * __restrict  __ATTR_ALIGN__(64) xre,
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
                        int32_t i;

                        for(i = 0; (i+95) < n; i += 96) {

                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                            zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0); 
                            zmm16 = _mm512_load_ps(&xre[i+64]);
                            zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm16,zmm17));
                            zmm18 = _mm512_load_ps(&xim[i+64]);
                            zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm18,zmm19));
                            zmm20 = _mm512_load_ps(&xre[i+80]);
                            zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm20,zmm21));
                            zmm22 = _mm512_load_ps(&xim[i+80]);
                            zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm22,zmm23));

                        }

                       for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));
                            zmm8  = _mm512_load_ps(&xre[i+32]);
                            zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm8,zmm9));
                            zmm10 = _mm512_load_ps(&xim[i+32]);
                             zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm10,zmm11));
                            zmm12 = _mm512_load_ps(&xre[i+48]);
                            zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm12,zmm13));
                            zmm14 = _mm512_load_ps(&xim[i+48]);
                            zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm14,zmm15));

                       }

                       for(; (i+31) < n; i += 32) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));
                            zmm4  = _mm512_load_ps(&xre[i+16]);
                            zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm4,zmm5));
                            zmm6  = _mm512_load_ps(&xim[i+16]);
                            zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm6,zmm7));

                       }

                       for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_load_ps(&xre[i+0]);
                            zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmm1));
                            zmm2  = _mm512_load_ps(&xim[i+0]);
                            zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm2,zmm3));

                       }

                        for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] + yre[i];
                            zim[i] = xim[i] + yim[i];
                      }

                  }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_16x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                          register __m512 zmm0,zmm1,zmm2,zmm3;
                          register __m512 zmm4,zmm5,zmm6,zmm7;
                          register __m512 zmm8,zmm9,zmm10,zmm11;
                          register __m512 zmm12,zmm13,zmm14,zmm15;
                          register __m512 zmm16,zmm17,zmm18,zmm19;
                          register __m512 zmm20,zmm21,zmm22,zmm23;
                          register __m512 zmm24,zmm25,zmm26,zmm27;
                          register __m512 zmm28,zmm29,zmm30,zmm31;
                          register __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+255) < n; i += 256) {
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                               zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                               zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                               zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                               zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                               zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                               zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                               zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                               zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                               zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                               zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                               zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                               zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                               zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                               zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_add_ps(zmm16,zmmx));
                               zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_add_ps(zmm17,zmmx));
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                               zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_add_ps(zmm18,zmmx));
                               zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_add_ps(zmm19,zmmx));
                               zmm20 = _mm512_loadu_ps(&xre[i+160]);
                              _mm512_storeu_ps(&zre[i+160], _mm512_add_ps(zmm20,zmmx));
                               zmm21 = _mm512_loadu_ps(&xim[i+160]);
                              _mm512_storeu_ps(&zim[i+160], _mm512_add_ps(zmm21,zmmx));
                               zmm22 = _mm512_loadu_ps(&xre[i+176]);
                              _mm512_storeu_ps(&zre[i+176], _mm512_add_ps(zmm22,zmmx));
                               zmm23 = _mm512_loadu_ps(&xim[i+176]);
                              _mm512_storeu_ps(&zim[i+176], _mm512_add_ps(zmm23,zmmx));
                               zmm24 = _mm512_loadu_ps(&xre[i+192]);
                              _mm512_storeu_ps(&zre[i+192], _mm512_add_ps(zmm24,zmmx));
                               zmm25 = _mm512_loadu_ps(&xim[i+192]);
                              _mm512_storeu_ps(&zim[i+192], _mm512_add_ps(zmm25,zmmx));
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                               zmm26 = _mm512_loadu_ps(&xre[i+208]);
                              _mm512_storeu_ps(&zre[i+208], _mm512_add_ps(zmm26,zmmx));
                               zmm27 = _mm512_loadu_ps(&xim[i+208]);
                              _mm512_storeu_ps(&zim[i+208], _mm512_add_ps(zmm27,zmmx));
                               zmm28 = _mm512_loadu_ps(&xre[i+224]);
                              _mm512_storeu_ps(&zre[i+224], _mm512_add_ps(zmm28,zmmx));
                               zmm29 = _mm512_loadu_ps(&xim[i+224]);
                              _mm512_storeu_ps(&zim[i+224], _mm512_add_ps(zmm29,zmmx));
                               zmm30 = _mm512_loadu_ps(&xre[i+240]);
                              _mm512_storeu_ps(&zre[i+240], _mm512_add_ps(zmm30,zmmx));
                               zmm31 = _mm512_loadu_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zim[i+240], _mm512_add_ps(zmm31,zmmx)); 

                         }

                        for(; (i+191) < n; i += 192) {
                               zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                               zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                               zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                               zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                               zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                               zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                               zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                               zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                               zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                               zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                               zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                               zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                               zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_add_ps(zmm16,zmmx));
                               zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_add_ps(zmm17,zmmx));
                               zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_add_ps(zmm18,zmmx));
                               zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_add_ps(zmm19,zmmx));
                               zmm20 = _mm512_loadu_ps(&xre[i+160]);
                              _mm512_storeu_ps(&zre[i+160], _mm512_add_ps(zmm20,zmmx));
                               zmm21 = _mm512_loadu_ps(&xim[i+160]);
                              _mm512_storeu_ps(&zim[i+160], _mm512_add_ps(zmm21,zmmx));
                               zmm22 = _mm512_loadu_ps(&xre[i+176]);
                              _mm512_storeu_ps(&zre[i+176], _mm512_add_ps(zmm22,zmmx));
                               zmm23 = _mm512_loadu_ps(&xim[i+176]);
                              _mm512_storeu_ps(&zim[i+176], _mm512_add_ps(zmm23,zmmx));

                        }

                        for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));

                        }

                       for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                             zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                             zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                             zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                             zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                             zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                             zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                             zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));

                       }

                      for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                             zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                             zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                             zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                      }

                     for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                             zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                     }
                    
                    
                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] + s;
                          zim[i] = zim[i] + s;
                    }
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                          register __m512 zmm0,zmm1,zmm2,zmm3;
                          register __m512 zmm4,zmm5,zmm6,zmm7;
                          register __m512 zmm8,zmm9,zmm10,zmm11;
                          register __m512 zmm12,zmm13,zmm14,zmm15;
                          register __m512 zmm16,zmm17,zmm18,zmm19;
                          register __m512 zmm20,zmm21,zmm22,zmm23;
                          register __m512 zmm24,zmm25,zmm26,zmm27;
                          register __m512 zmm28,zmm29,zmm30,zmm31;
                          register __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+255) < n; i += 256) {
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                              zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_add_ps(zmm16,zmmx));
                              zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_add_ps(zmm17,zmmx));
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                              zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_add_ps(zmm18,zmmx));
                              zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_add_ps(zmm19,zmmx));
                              zmm20 = _mm512_load_ps(&xre[i+160]);
                              _mm512_store_ps(&zre[i+160], _mm512_add_ps(zmm20,zmmx));
                              zmm21 = _mm512_load_ps(&xim[i+160]);
                              _mm512_store_ps(&zim[i+160], _mm512_add_ps(zmm21,zmmx));
                              zmm22 = _mm512_load_ps(&xre[i+176]);
                              _mm512_store_ps(&zre[i+176], _mm512_add_ps(zmm22,zmmx));
                              zmm23 = _mm512_load_ps(&xim[i+176]);
                              _mm512_store_ps(&zim[i+176], _mm512_add_ps(zmm23,zmmx));
                              zmm24 = _mm512_load_ps(&xre[i+192]);
                              _mm512_store_ps(&zre[i+192], _mm512_add_ps(zmm24,zmmx));
                              zmm25 = _mm512_load_ps(&xim[i+192]);
                              _mm512_storeu_ps(&zim[i+192], _mm512_add_ps(zmm25,zmmx));
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                              zmm26 = _mm512_load_ps(&xre[i+208]);
                              _mm512_store_ps(&zre[i+208], _mm512_add_ps(zmm26,zmmx));
                              zmm27 = _mm512_load_ps(&xim[i+208]);
                              _mm512_store_ps(&zim[i+208], _mm512_add_ps(zmm27,zmmx));
                              zmm28 = _mm512_load_ps(&xre[i+224]);
                              _mm512_store_ps(&zre[i+224], _mm512_add_ps(zmm28,zmmx));
                              zmm29 = _mm512_load_ps(&xim[i+224]);
                              _mm512_store_ps(&zim[i+224], _mm512_add_ps(zmm29,zmmx));
                              zmm30 = _mm512_load_ps(&xre[i+240]);
                              _mm512_store_ps(&zre[i+240], _mm512_add_ps(zmm30,zmmx));
                              zmm31 = _mm512_load_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zim[i+240], _mm512_add_ps(zmm31,zmmx)); 

                         }

                        for(; (i+191) < n; i += 192) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                              zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_add_ps(zmm16,zmmx));
                              zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_add_ps(zmm17,zmmx));
                              zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_add_ps(zmm18,zmmx));
                              zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_add_ps(zmm19,zmmx));
                              zmm20 = _mm512_load_ps(&xre[i+160]);
                              _mm512_store_ps(&zre[i+160], _mm512_add_ps(zmm20,zmmx));
                              zmm21 = _mm512_load_ps(&xim[i+160]);
                              _mm512_store_ps(&zim[i+160], _mm512_add_ps(zmm21,zmmx));
                              zmm22 = _mm512_load_ps(&xre[i+176]);
                              _mm512_store_ps(&zre[i+176], _mm512_add_ps(zmm22,zmmx));
                              zmm23 = _mm512_load_ps(&xim[i+176]);
                              _mm512_store_ps(&zim[i+176], _mm512_add_ps(zmm23,zmmx)); 

                        }

                        for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                              zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                        }

                       for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                       }

                      for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                      }

                     for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                     }
                    
                    
                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] + s;
                          zim[i] = zim[i] + s;
                    }
                 }




                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                          register __m512 zmm0,zmm1,zmm2,zmm3;
                          register __m512 zmm4,zmm5,zmm6,zmm7;
                          register __m512 zmm8,zmm9,zmm10,zmm11;
                          register __m512 zmm12,zmm13,zmm14,zmm15;
                          register __m512 zmm16,zmm17,zmm18,zmm19;
                          register  __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                              zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_add_ps(zmm16,zmmx));
                              zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_add_ps(zmm17,zmmx));
                              zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_add_ps(zmm18,zmmx));
                              zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_add_ps(zmm19,zmmx));
                        }

                       for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                              zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                       }

                       for(; (i+79) < n; i += 80) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                       }

                      for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx)); 
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));    
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] + s;
                          zim[i] = zim[i] + s;
                    }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                               vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                          register __m512 zmm0,zmm1,zmm2,zmm3;
                          register __m512 zmm4,zmm5,zmm6,zmm7;
                          register __m512 zmm8,zmm9,zmm10,zmm11;
                          register __m512 zmm12,zmm13,zmm14,zmm15;
                          register __m512 zmm16,zmm17,zmm18,zmm19;
                          register  zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                               zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                               zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                               zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_add_ps(zmm16,zmmx));
                               zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_add_ps(zmm17,zmmx));
                               zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_add_ps(zmm18,zmmx));
                               zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_add_ps(zmm19,zmmx));
                        }

                       for(; (i+127) < n; i += 128) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                               zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_add_ps(zmm12,zmmx));
                              zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_add_ps(zmm13,zmmx));
                              zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_add_ps(zmm14,zmmx));
                              zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_add_ps(zmm15,zmmx));
                       }

                       for(; (i+79) < n; i += 80) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                       }

                      for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                               zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] + s;
                          zim[i] = zim[i] + s;
                    }

                 }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                          register __m512 zmm0,zmm1,zmm2,zmm3;
                          register __m512 zmm4,zmm5,zmm6,zmm7;
                          register __m512 zmm8,zmm9,zmm10,zmm11;
                          register  __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                              zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                              zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                              zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                            
                        }

                       for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                               zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                      }

                     for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                               zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                              zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));  
                     }

                    for(; (i+15) < n; i += 16) {
                              zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                               zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] + s;
                          zim[i] = zim[i] + s;
                    }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cadd_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                          register __m512 zmm0,zmm1,zmm2,zmm3;
                          register __m512 zmm4,zmm5,zmm6,zmm7;
                          register __m512 zmm8,zmm9,zmm10,zmm11;
                          register  __m512 zmmx    = vs;
                          const float * __restrict pzmm =  (const float*)&vs[0];
                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                               zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                               zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_add_ps(zmm8,zmmx));
                               zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_add_ps(zmm9,zmmx));
                              zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_add_ps(zmm10,zmmx));
                               zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_add_ps(zmm11,zmmx));
                        }

                        for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                              zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                              zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                               zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                              zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_add_ps(zmm4,zmmx));
                              zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_add_ps(zmm5,zmmx));
                              zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_add_ps(zmm6,zmmx));
                              zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_add_ps(zmm7,zmmx));
                      }

                     for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                             zmm1 = _mm512_load_ps(&xim[i+0]);
                             _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                             zmm2 = _mm512_load_ps(&xre[i+16]);
                             _mm512_store_ps(&zre[i+16], _mm512_add_ps(zmm2,zmmx));
                             zmm3 = _mm512_load_ps(&xim[i+16]);
                             _mm512_store_ps(&zim[i+16], _mm512_add_ps(zmm3,zmmx));
                     }

                    for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             _mm512_store_ps(&zre[i+0], _mm512_add_ps(zmm0,zmmx));
                             zmm1 = _mm512_load_ps(&xim[i+0]);
                             _mm512_store_ps(&zim[i+0], _mm512_add_ps(zmm1,zmmx));
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] + s;
                          zim[i] = zim[i] + s;
                    }

                 }


                  

        } // math

} // gms















#endif /*__GMS_CADD_VEC_ZMM16R4_HPP__*/
