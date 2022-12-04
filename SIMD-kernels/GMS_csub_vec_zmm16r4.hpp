

#ifndef __GMS_CSUB_VEC_ZMM16R4_HPP__
#define __GMS_CSUB_VEC_ZMM16R4_HPP__


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

    const unsigned int GMS_CSUB_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CSUB_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CSUB_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CSUB_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CSUB_VEC_ZMM16R4_MAJOR+
      100U*GMS_CSUB_VEC_ZMM16R4_MINOR+
      10U*GMS_CSUB_VEC_ZMM16R4_MICRO;
    const char * const GMS_CSUB_VEC_ZMM16R4_CREATION_DATE = "01-12-2022 15:59 AM +00200 (THR 01 DEC 2022 GMT+2)";
    const char * const GMS_CSUB_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CSUB_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CSUB_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector subtraction operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace  gms {

           namespace  math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_16x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            const __m512 zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            const __m512 zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            const __m512 zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            const __m512 zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            const __m512 zmm32 = _mm512_loadu_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            const __m512 zmm34 = _mm512_loadu_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            const __m512 zmm36 = _mm512_loadu_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            const __m512 zmm38 = _mm512_loadu_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
                            const __m512 zmm40 = _mm512_loadu_ps(&xre[i+160]); //SPILLING!!
                            const __m512 zmm41 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm40,zmm41));
                            const __m512 zmm42 = _mm512_loadu_ps(&xim[i+160]);
                            const __m512 zmm43 = _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&zim[i+160], _mm512_sub_ps(zmm42,zmm43));
                            const __m512 zmm44 = _mm512_loadu_ps(&xre[i+176]); //SPILLING!!
                            const __m512 zmm45 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm44,zmm45));
                            const __m512 zmm46 = _mm512_loadu_ps(&xim[i+176]);
                            const __m512 zmm47 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&zim[i+176], _mm512_sub_ps(zmm46,zmm47));
                            const __m512 zmm48 = _mm512_loadu_ps(&xre[i+192]); //SPILLING!!
                            const __m512 zmm49 = _mm512_loadu_ps(&yre[i+192]);
                            _mm512_storeu_ps(&zre[i+192], _mm512_sub_ps(zmm48,zmm49));
                            const __m512 zmm50 = _mm512_loadu_ps(&xim[i+192]);
                            const __m512 zmm51 = _mm512_loadu_ps(&yim[i+192]);
                            _mm512_storeu_ps(&zim[i+192], _mm512_sub_ps(zmm50,zmm51));
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+255],_MM_HINT_T0);
                            const __m512 zmm52 = _mm512_loadu_ps(&xre[i+208]); //SPILLING!!
                            const __m512 zmm53 = _mm512_loadu_ps(&yre[i+208]);
                            _mm512_storeu_ps(&zre[i+208], _mm512_sub_ps(zmm52,zmm53));
                            const __m512 zmm54 = _mm512_loadu_ps(&xim[i+208]);
                            const __m512 zmm55 = _mm512_loadu_ps(&yim[i+208]);
                            _mm512_storeu_ps(&zim[i+208], _mm512_sub_ps(zmm54,zmm55));
                            const __m512 zmm56 = _mm512_loadu_ps(&xre[i+224]); //SPILLING!!
                            const __m512 zmm57 = _mm512_loadu_ps(&yre[i+224]);
                            _mm512_storeu_ps(&zre[i+224], _mm512_sub_ps(zmm56,zmm57));
                            const __m512 zmm58 = _mm512_loadu_ps(&xim[i+224]);
                            const __m512 zmm59 = _mm512_loadu_ps(&yim[i+224]);
                            _mm512_storeu_ps(&zim[i+224], _mm512_sub_ps(zmm58,zmm59));
                            const __m512 zmm60 = _mm512_loadu_ps(&xre[i+240]); //SPILLING!!
                            const __m512 zmm61 = _mm512_loadu_ps(&yre[i+240]);
                            _mm512_storeu_ps(&zre[i+240], _mm512_sub_ps(zmm60,zmm61));
                            const __m512 zmm62 = _mm512_loadu_ps(&xim[i+240]);
                            const __m512 zmm63 = _mm512_loadu_ps(&yim[i+240]);
                            _mm512_storeu_ps(&zim[i+240], _mm512_sub_ps(zmm62,zmm63));
                            
#else
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+255],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            const __m512 zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            const __m512 zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            const __m512 zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            const __m512 zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            const __m512 zmm32 = _mm512_loadu_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_loadu_ps(&yre[i+128]);
                            const __m512 zmm34 = _mm512_loadu_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_loadu_ps(&yim[i+128]);
                            const __m512 zmm36 = _mm512_loadu_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_loadu_ps(&yre[i+144]);
                            const __m512 zmm38 = _mm512_loadu_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_loadu_ps(&yim[i+144]);
                            const __m512 zmm40 = _mm512_loadu_ps(&xre[i+160]); //SPILLING!!
                            const __m512 zmm41 = _mm512_loadu_ps(&yre[i+160]);
                            const __m512 zmm42 = _mm512_loadu_ps(&xim[i+160]);
                            const __m512 zmm43 = _mm512_loadu_ps(&yim[i+160]);
                            const __m512 zmm44 = _mm512_loadu_ps(&xre[i+176]); //SPILLING!!
                            const __m512 zmm45 = _mm512_loadu_ps(&yre[i+176]);
                            const __m512 zmm46 = _mm512_loadu_ps(&xim[i+176]);
                            const __m512 zmm47 = _mm512_loadu_ps(&yim[i+176]);
                            const __m512 zmm48 = _mm512_loadu_ps(&xre[i+192]); //SPILLING!!
                            const __m512 zmm49 = _mm512_loadu_ps(&yre[i+192]);
                            const __m512 zmm50 = _mm512_loadu_ps(&xim[i+192]);
                            const __m512 zmm51 = _mm512_loadu_ps(&yim[i+192]);
                            const __m512 zmm52 = _mm512_loadu_ps(&xre[i+208]); //SPILLING!!
                            const __m512 zmm53 = _mm512_loadu_ps(&yre[i+208]);
                            const __m512 zmm54 = _mm512_loadu_ps(&xim[i+208]);
                            const __m512 zmm55 = _mm512_loadu_ps(&yim[i+208]);
                            const __m512 zmm56 = _mm512_loadu_ps(&xre[i+224]); //SPILLING!!
                            const __m512 zmm57 = _mm512_loadu_ps(&yre[i+224]);
                            const __m512 zmm58 = _mm512_loadu_ps(&xim[i+224]);
                            const __m512 zmm59 = _mm512_loadu_ps(&yim[i+224]);
                            const __m512 zmm60 = _mm512_loadu_ps(&xre[i+240]); //SPILLING!!
                            const __m512 zmm61 = _mm512_loadu_ps(&yre[i+240]);
                            const __m512 zmm62 = _mm512_loadu_ps(&xim[i+240]);
                            const __m512 zmm63 = _mm512_loadu_ps(&yim[i+240]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
                            _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm40,zmm41));
                            _mm512_storeu_ps(&zim[i+160], _mm512_sub_ps(zmm42,zmm43));
                            _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm44,zmm45));
                            _mm512_storeu_ps(&zim[i+176], _mm512_sub_ps(zmm46,zmm47));
                            _mm512_storeu_ps(&zre[i+192], _mm512_sub_ps(zmm48,zmm49));
                            _mm512_storeu_ps(&zim[i+192], _mm512_sub_ps(zmm50,zmm51));
                            _mm512_storeu_ps(&zre[i+208], _mm512_sub_ps(zmm52,zmm53));
                            _mm512_storeu_ps(&zim[i+208], _mm512_sub_ps(zmm54,zmm55));
                            _mm512_storeu_ps(&zre[i+224], _mm512_sub_ps(zmm56,zmm57));
                            _mm512_storeu_ps(&zim[i+224], _mm512_sub_ps(zmm58,zmm59));
                            _mm512_storeu_ps(&zre[i+240], _mm512_sub_ps(zmm60,zmm61));
                            _mm512_storeu_ps(&zim[i+240], _mm512_sub_ps(zmm62,zmm63));
#endif
                        }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1                      
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            const __m512 zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            const __m512 zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            const __m512 zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            const __m512 zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            const __m512 zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            const __m512 zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            const __m512 zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            const __m512 zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1                      
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1    
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#endif
                       }

                      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1    
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
#endif
                      }

                      for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_16x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yim,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                        int32_t i;
                        
                        for(i = 0; (i+255) < n; i += 256) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            const __m512 zmm24 = _mm512_load_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            const __m512 zmm26 = _mm512_load_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            const __m512 zmm28 = _mm512_load_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            const __m512 zmm30 = _mm512_load_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            const __m512 zmm32 = _mm512_load_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            const __m512 zmm34 = _mm512_load_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            const __m512 zmm36 = _mm512_load_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            const __m512 zmm38 = _mm512_load_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
                            const __m512 zmm40 = _mm512_load_ps(&xre[i+160]); //SPILLING!!
                            const __m512 zmm41 = _mm512_load_ps(&yre[i+160]);
                            _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm40,zmm41));
                            const __m512 zmm42 = _mm512_load_ps(&xim[i+160]);
                            const __m512 zmm43 = _mm512_load_ps(&yim[i+160]);
                            _mm512_store_ps(&zim[i+160], _mm512_sub_ps(zmm42,zmm43));
                            const __m512 zmm44 = _mm512_load_ps(&xre[i+176]); //SPILLING!!
                            const __m512 zmm45 = _mm512_load_ps(&yre[i+176]);
                            _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm44,zmm45));
                            const __m512 zmm46 = _mm512_load_ps(&xim[i+176]);
                            const __m512 zmm47 = _mm512_load_ps(&yim[i+176]);
                            _mm512_store_ps(&zim[i+176], _mm512_sub_ps(zmm46,zmm47));
                            const __m512 zmm48 = _mm512_load_ps(&xre[i+192]); //SPILLING!!
                            const __m512 zmm49 = _mm512_load_ps(&yre[i+192]);
                            _mm512_store_ps(&zre[i+192], _mm512_sub_ps(zmm48,zmm49));
                            const __m512 zmm50 = _mm512_load_ps(&xim[i+192]);
                            const __m512 zmm51 = _mm512_load_ps(&yim[i+192]);
                            _mm512_store_ps(&zim[i+192], _mm512_sub_ps(zmm50,zmm51));
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+255],_MM_HINT_T0);
                            const __m512 zmm52 = _mm512_load_ps(&xre[i+208]); //SPILLING!!
                            const __m512 zmm53 = _mm512_load_ps(&yre[i+208]);
                            _mm512_store_ps(&zre[i+208], _mm512_sub_ps(zmm52,zmm53));
                            const __m512 zmm54 = _mm512_load_ps(&xim[i+208]);
                            const __m512 zmm55 = _mm512_load_ps(&yim[i+208]);
                            _mm512_store_ps(&zim[i+208], _mm512_sub_ps(zmm54,zmm55));
                            const __m512 zmm56 = _mm512_load_ps(&xre[i+224]); //SPILLING!!
                            const __m512 zmm57 = _mm512_load_ps(&yre[i+224]);
                            _mm512_store_ps(&zre[i+224], _mm512_sub_ps(zmm56,zmm57));
                            const __m512 zmm58 = _mm512_load_ps(&xim[i+224]);
                            const __m512 zmm59 = _mm512_load_ps(&yim[i+224]);
                            _mm512_store_ps(&zim[i+224], _mm512_sub_ps(zmm58,zmm59));
                            const __m512 zmm60 = _mm512_load_ps(&xre[i+240]); //SPILLING!!
                            const __m512 zmm61 = _mm512_load_ps(&yre[i+240]);
                            _mm512_store_ps(&zre[i+240], _mm512_sub_ps(zmm60,zmm61));
                            const __m512 zmm62 = _mm512_load_ps(&xim[i+240]);
                            const __m512 zmm63 = _mm512_load_ps(&yim[i+240]);
                            _mm512_store_ps(&zim[i+240], _mm512_sub_ps(zmm62,zmm63));
                            
#else
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+255],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            const __m512 zmm24 = _mm512_load_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_load_ps(&yre[i+96]);
                            const __m512 zmm26 = _mm512_load_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_load_ps(&yim[i+96]);
                            const __m512 zmm28 = _mm512_load_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_load_ps(&yre[i+112]);
                            const __m512 zmm30 = _mm512_load_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_load_ps(&yim[i+112]);
                            const __m512 zmm32 = _mm512_load_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_load_ps(&yre[i+128]);
                            const __m512 zmm34 = _mm512_load_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_load_ps(&yim[i+128]);
                            const __m512 zmm36 = _mm512_load_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_load_ps(&yre[i+144]);
                            const __m512 zmm38 = _mm512_load_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_load_ps(&yim[i+144]);
                            const __m512 zmm40 = _mm512_load_ps(&xre[i+160]); //SPILLING!!
                            const __m512 zmm41 = _mm512_load_ps(&yre[i+160]);
                            const __m512 zmm42 = _mm512_load_ps(&xim[i+160]);
                            const __m512 zmm43 = _mm512_load_ps(&yim[i+160]);
                            const __m512 zmm44 = _mm512_load_ps(&xre[i+176]); //SPILLING!!
                            const __m512 zmm45 = _mm512_load_ps(&yre[i+176]);
                            const __m512 zmm46 = _mm512_load_ps(&xim[i+176]);
                            const __m512 zmm47 = _mm512_load_ps(&yim[i+176]);
                            const __m512 zmm48 = _mm512_load_ps(&xre[i+192]); //SPILLING!!
                            const __m512 zmm49 = _mm512_load_ps(&yre[i+192]);
                            const __m512 zmm50 = _mm512_load_ps(&xim[i+192]);
                            const __m512 zmm51 = _mm512_load_ps(&yim[i+192]);
                            const __m512 zmm52 = _mm512_load_ps(&xre[i+208]); //SPILLING!!
                            const __m512 zmm53 = _mm512_load_ps(&yre[i+208]);
                            const __m512 zmm54 = _mm512_load_ps(&xim[i+208]);
                            const __m512 zmm55 = _mm512_load_ps(&yim[i+208]);
                            const __m512 zmm56 = _mm512_load_ps(&xre[i+224]); //SPILLING!!
                            const __m512 zmm57 = _mm512_load_ps(&yre[i+224]);
                            const __m512 zmm58 = _mm512_load_ps(&xim[i+224]);
                            const __m512 zmm59 = _mm512_load_ps(&yim[i+224]);
                            const __m512 zmm60 = _mm512_load_ps(&xre[i+240]); //SPILLING!!
                            const __m512 zmm61 = _mm512_load_ps(&yre[i+240]);
                            const __m512 zmm62 = _mm512_load_ps(&xim[i+240]);
                            const __m512 zmm63 = _mm512_load_ps(&yim[i+240]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
                            _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm40,zmm41));
                            _mm512_store_ps(&zim[i+160], _mm512_sub_ps(zmm42,zmm43));
                            _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm44,zmm45));
                            _mm512_store_ps(&zim[i+176], _mm512_sub_ps(zmm46,zmm47));
                            _mm512_store_ps(&zre[i+192], _mm512_sub_ps(zmm48,zmm49));
                            _mm512_store_ps(&zim[i+192], _mm512_sub_ps(zmm50,zmm51));
                            _mm512_store_ps(&zre[i+208], _mm512_sub_ps(zmm52,zmm53));
                            _mm512_store_ps(&zim[i+208], _mm512_sub_ps(zmm54,zmm55));
                            _mm512_store_ps(&zre[i+224], _mm512_sub_ps(zmm56,zmm57));
                            _mm512_store_ps(&zim[i+224], _mm512_sub_ps(zmm58,zmm59));
                            _mm512_store_ps(&zre[i+240], _mm512_sub_ps(zmm60,zmm61));
                            _mm512_store_ps(&zim[i+240], _mm512_sub_ps(zmm62,zmm63));
#endif
                        }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1                      
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            const __m512 zmm24 = _mm512_load_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            const __m512 zmm26 = _mm512_load_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            const __m512 zmm28 = _mm512_load_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            const __m512 zmm30 = _mm512_load_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            const __m512 zmm24 = _mm512_load_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_load_ps(&yre[i+96]);
                            const __m512 zmm26 = _mm512_load_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_load_ps(&yim[i+96]);
                            const __m512 zmm28 = _mm512_load_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_load_ps(&yre[i+112]);
                            const __m512 zmm30 = _mm512_load_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1                      
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1    
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#endif
                       }

                      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1    
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
#endif
                      }

                      for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_10x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                        int32_t i;

                        for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            const __m512 zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            const __m512 zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            const __m512 zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            const __m512 zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            const __m512 zmm32 = _mm512_loadu_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            const __m512 zmm34 = _mm512_loadu_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            const __m512 zmm36 = _mm512_loadu_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            const __m512 zmm38 = _mm512_loadu_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
#else
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            const __m512 zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            const __m512 zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            const __m512 zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            const __m512 zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            const __m512 zmm32 = _mm512_loadu_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_loadu_ps(&yre[i+128]);
                            const __m512 zmm34 = _mm512_loadu_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_loadu_ps(&yim[i+128]);
                            const __m512 zmm36 = _mm512_loadu_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_loadu_ps(&yre[i+144]);
                            const __m512 zmm38 = _mm512_loadu_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
#endif
                        }

                       for(; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#endif
                      }

                      for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
#endif
                     }

                     for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_10x_a(const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yim,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                        int32_t i;

                        for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            const __m512 zmm24 = _mm512_load_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_load_ps(&yre[i+96]);
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            const __m512 zmm26 = _mm512_load_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_load_ps(&yim[i+96]);
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            const __m512 zmm28 = _mm512_load_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_load_ps(&yre[i+112]);
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            const __m512 zmm30 = _mm512_load_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_load_ps(&yim[i+112]);
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            const __m512 zmm32 = _mm512_load_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_load_ps(&yre[i+128]);
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            const __m512 zmm34 = _mm512_load_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_load_ps(&yim[i+128]);
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            const __m512 zmm36 = _mm512_load_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_load_ps(&yre[i+144]);
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            const __m512 zmm38 = _mm512_load_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
#else
                            _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+128],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            const __m512 zmm24 = _mm512_load_ps(&xre[i+96]);
                            const __m512 zmm25 = _mm512_load_ps(&yre[i+96]);
                            const __m512 zmm26 = _mm512_load_ps(&xim[i+96]);
                            const __m512 zmm27 = _mm512_load_ps(&yim[i+96]);
                            const __m512 zmm28 = _mm512_load_ps(&xre[i+112]);
                            const __m512 zmm29 = _mm512_load_ps(&yre[i+112]);
                            const __m512 zmm30 = _mm512_load_ps(&xim[i+112]);
                            const __m512 zmm31 = _mm512_load_ps(&yim[i+112]);
                            const __m512 zmm32 = _mm512_load_ps(&xre[i+128]); //SPILLING!!
                            const __m512 zmm33 = _mm512_load_ps(&yre[i+128]);
                            const __m512 zmm34 = _mm512_load_ps(&xim[i+128]);
                            const __m512 zmm35 = _mm512_load_ps(&yim[i+128]);
                            const __m512 zmm36 = _mm512_load_ps(&xre[i+144]); //SPILLING!!
                            const __m512 zmm37 = _mm512_load_ps(&yre[i+144]);
                            const __m512 zmm38 = _mm512_load_ps(&xim[i+144]);
                            const __m512 zmm39 = _mm512_load_ps(&yim[i+144]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
                            _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm24,zmm25));
                            _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm26,zmm27));
                            _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm28,zmm29));
                            _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm30,zmm31));
                            _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm32,zmm33));
                            _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm34,zmm35));
                            _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm36,zmm37));
                            _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm38,zmm39));
#endif
                        }

                       for(; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeups(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#endif
                      }

                      for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
#endif
                     }

                     for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_6x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                        int32_t i;

                        for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0); 
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#else
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0); 
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            const __m512 zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            const __m512 zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            const __m512 zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            const __m512 zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            const __m512 zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            const __m512 zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            const __m512 zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            const __m512 zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#endif
                        }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            const __m512 zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&yre[i+32]);
                            const __m512 zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            const __m512 zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            const __m512 zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            const __m512 zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            const __m512 zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            const __m512 zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_loadu_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#endif
                       }

                       for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_loadu_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_storeu_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
#endif
                       }

                        for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_6x_a( const float * __restrict  __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict  __ATTR_ALIGN__(64) yim,
                                                  float * __restrict        __ATTR_ALIGN__(64) zre,
                                                  float * __restrict        __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                        int32_t i;

                        for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            register const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            register const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            register const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            register const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            register const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            register const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            register const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            register const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            register const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            register const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            register const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            register const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0); 
                            register const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            register const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            register const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            register const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            register const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            register const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            register const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            register const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#else
                            _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0); 
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            register const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            register const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            register const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            register const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            register const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            register const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            register const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            register const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            register const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            register const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            register const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            register const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            register const __m512 zmm16 = _mm512_load_ps(&xre[i+64]);
                            register const __m512 zmm17 = _mm512_load_ps(&yre[i+64]);
                            register const __m512 zmm18 = _mm512_load_ps(&xim[i+64]);
                            register const __m512 zmm19 = _mm512_load_ps(&yim[i+64]);
                            register const __m512 zmm20 = _mm512_load_ps(&xre[i+80]);
                            register const __m512 zmm21 = _mm512_load_ps(&yre[i+80]);
                            register const __m512 zmm22 = _mm512_load_ps(&xim[i+80]);
                            register const __m512 zmm23 = _mm512_load_ps(&yim[i+80]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
                            _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm16,zmm17));
                            _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm18,zmm19));
                            _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm20,zmm21));
                            _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm22,zmm23));
#endif
                        }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            register const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            register const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            register const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            register const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            register const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            register const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            register const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            register const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            register const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            register const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            register const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            register const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#else
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            register const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            register const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            register const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            register const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            register const __m512 zmm8  = _mm512_load_ps(&xre[i+32]);
                            register const __m512 zmm9  = _mm512_load_ps(&yre[i+32]);
                            register const __m512 zmm10 = _mm512_load_ps(&xim[i+32]);
                            register const __m512 zmm11 = _mm512_load_ps(&yim[i+32]);
                            register const __m512 zmm12 = _mm512_load_ps(&xre[i+48]);
                            register const __m512 zmm13 = _mm512_load_ps(&yre[i+48]);
                            register const __m512 zmm14 = _mm512_load_ps(&xim[i+48]);
                            register const __m512 zmm15 = _mm512_load_ps(&yim[i+48]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
                            _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm8,zmm9));
                            _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm10,zmm11));
                            _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm12,zmm13));
                            _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm14,zmm15));
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
                            register const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            register const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            register const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            register const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#else
                            const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            const __m512 zmm4  = _mm512_load_ps(&xre[i+16]);
                            const __m512 zmm5  = _mm512_load_ps(&yre[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&xim[i+16]);
                            const __m512 zmm7  = _mm512_load_ps(&yim[i+16]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
                            _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm4,zmm5));
                            _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm6,zmm7));
#endif
                       }

                       for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1   
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmm1));
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm2,zmm3));
#else
                            register const __m512 zmm0  = _mm512_load_ps(&xre[i+0]);
                            register const __m512 zmm1  = _mm512_load_ps(&yre[i+0]);
                            register const __m512 zmm2  = _mm512_load_ps(&xim[i+0]);
                            register const __m512 zmm3  = _mm512_load_ps(&yim[i+0]);
                            _mm512_store_ps(&zre[i+0],  _mm512_sub_ps(zmm0,zmm1));
                            _mm512_store_ps(&zim[i+0],  _mm512_sub_ps(zmm2,zmm3));
#endif
                       }

                        for(; (i+0) < n; i += 1) {
                             
                            zre[i] = xre[i] - yre[i];
                            zim[i] = xim[i] - yim[i];
                      }

                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_16x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                         register const __m512 zmmx    = vs;
                         const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+255) < n; i += 256) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              register const __m512 zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              register const __m512 zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              register const __m512 zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              register const __m512 zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              register const __m512 zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              register const __m512 zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                              register const __m512 zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              register const __m512 zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
                              register const __m512 zmm20 = _mm512_loadu_ps(&xre[i+160]);
                              _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm20,zmmx));
                              register const __m512 zmm21 = _mm512_loadu_ps(&xim[i+160]);
                              _mm512_storeu_ps(&zim[i+160], _mm512_sub_ps(zmm21,zmmx));
                              register const __m512 zmm22 = _mm512_loadu_ps(&xre[i+176]);
                              _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm22,zmmx));
                              register const __m512 zmm23 = _mm512_loadu_ps(&xim[i+176]);
                              _mm512_storeu_ps(&zim[i+176], _mm512_sub_ps(zmm23,zmmx));
                              register const __m512 zmm24 = _mm512_loadu_ps(&xre[i+192]);
                              _mm512_storeu_ps(&zre[i+192], _mm512_sub_ps(zmm24,zmmx));
                              register const __m512 zmm25 = _mm512_loadu_ps(&xim[i+192]);
                              _mm512_storeu_ps(&zim[i+192], _mm512_sub_ps(zmm25,zmmx));
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                              register const __m512 zmm26 = _mm512_loadu_ps(&xre[i+208]);
                              _mm512_storeu_ps(&zre[i+208], _mm512_sub_ps(zmm26,zmmx));
                              register const __m512 zmm27 = _mm512_loadu_ps(&xim[i+208]);
                              _mm512_storeu_ps(&zim[i+208], _mm512_sub_ps(zmm27,zmmx));
                              register const __m512 zmm28 = _mm512_loadu_ps(&xre[i+224]);
                              _mm512_storeu_ps(&zre[i+224], _mm512_sub_ps(zmm28,zmmx));
                              register const __m512 zmm29 = _mm512_loadu_ps(&xim[i+224]);
                              _mm512_storeu_ps(&zim[i+224], _mm512_sub_ps(zmm29,zmmx));
                              register const __m512 zmm30 = _mm512_loadu_ps(&xre[i+240]);
                              _mm512_storeu_ps(&zre[i+240], _mm512_sub_ps(zmm30,zmmx));
                              register const __m512 zmm31 = _mm512_loadu_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zim[i+240], _mm512_sub_ps(zmm31,zmmx)); 
#else
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              register const __m512 zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              register const __m512 zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              register const __m512 zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              register const __m512 zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              register const __m512 zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              register const __m512 zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              register const __m512 zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              register const __m512 zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              register const __m512 zmm20 = _mm512_loadu_ps(&xre[i+160]);
                              register const __m512 zmm21 = _mm512_loadu_ps(&xim[i+160]);
                              register const __m512 zmm22 = _mm512_loadu_ps(&xre[i+176]);
                              register const __m512 zmm23 = _mm512_loadu_ps(&xim[i+176]);
                              register const __m512 zmm24 = _mm512_loadu_ps(&xre[i+192]);
                              register const __m512 zmm25 = _mm512_loadu_ps(&xim[i+192]);
                              register const __m512 zmm26 = _mm512_loadu_ps(&xre[i+208]);
                              register const __m512 zmm27 = _mm512_loadu_ps(&xim[i+208]);
                              register const __m512 zmm28 = _mm512_loadu_ps(&xre[i+224]);
                              register const __m512 zmm29 = _mm512_loadu_ps(&xim[i+224]);
                              register const __m512 zmm30 = _mm512_loadu_ps(&xre[i+240]);
                              register const __m512 zmm31 = _mm512_loadu_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
                              _mm512_storeu_ps(&zre[i+160], _mm512_sub_ps(zmm20,zmmx));
                              _mm512_storeu_ps(&zim[i+160], _mm512_sub_ps(zmm21,zmmx));
                              _mm512_storeu_ps(&zre[i+176], _mm512_sub_ps(zmm22,zmmx));
                              _mm512_storeu_ps(&zim[i+176], _mm512_sub_ps(zmm23,zmmx));
                              _mm512_storeu_ps(&zre[i+192], _mm512_sub_ps(zmm24,zmmx));
                              _mm512_storeu_ps(&zim[i+192], _mm512_sub_ps(zmm25,zmmx));
                              _mm512_storeu_ps(&zre[i+208], _mm512_sub_ps(zmm26,zmmx));
                              _mm512_storeu_ps(&zim[i+208], _mm512_sub_ps(zmm27,zmmx));
                              _mm512_storeu_ps(&zre[i+224], _mm512_sub_ps(zmm28,zmmx));
                              _mm512_storeu_ps(&zim[i+224], _mm512_sub_ps(zmm29,zmmx));
                              _mm512_storeu_ps(&zre[i+240], _mm512_sub_ps(zmm30,zmmx));
                              _mm512_storeu_ps(&zim[i+240], _mm512_sub_ps(zmm31,zmmx)); 
                     
#endif
                         }

                        for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              register const __m512 zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              register const __m512 zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              register const __m512 zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              register const __m512 zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              register const __m512 zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              register const __m512 zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              register const __m512 zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              register const __m512 zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
#endif
                        }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#endif
                       }

                      for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#endif
                     }
                    
                    
                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] - s;
                          zim[i] = zim[i] - s;
                    }
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                         register const __m512 zmmx    = vs;
                         const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+255) < n; i += 256) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              register const __m512 zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              register const __m512 zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              register const __m512 zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              register const __m512 zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              register const __m512 zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              register const __m512 zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                              register const __m512 zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              register const __m512 zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
                              register const __m512 zmm20 = _mm512_load_ps(&xre[i+160]);
                              _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm20,zmmx));
                              register const __m512 zmm21 = _mm512_load_ps(&xim[i+160]);
                              _mm512_store_ps(&zim[i+160], _mm512_sub_ps(zmm21,zmmx));
                              register const __m512 zmm22 = _mm512_load_ps(&xre[i+176]);
                              _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm22,zmmx));
                              register const __m512 zmm23 = _mm512_load_ps(&xim[i+176]);
                              _mm512_store_ps(&zim[i+176], _mm512_sub_ps(zmm23,zmmx));
                              register const __m512 zmm24 = _mm512_load_ps(&xre[i+192]);
                              _mm512_store_ps(&zre[i+192], _mm512_sub_ps(zmm24,zmmx));
                              register const __m512 zmm25 = _mm512_load_ps(&xim[i+192]);
                              _mm512_storeu_ps(&zim[i+192], _mm512_sub_ps(zmm25,zmmx));
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                              register const __m512 zmm26 = _mm512_load_ps(&xre[i+208]);
                              _mm512_store_ps(&zre[i+208], _mm512_sub_ps(zmm26,zmmx));
                              register const __m512 zmm27 = _mm512_load_ps(&xim[i+208]);
                              _mm512_store_ps(&zim[i+208], _mm512_sub_ps(zmm27,zmmx));
                              register const __m512 zmm28 = _mm512_load_ps(&xre[i+224]);
                              _mm512_store_ps(&zre[i+224], _mm512_sub_ps(zmm28,zmmx));
                              register const __m512 zmm29 = _mm512_load_ps(&xim[i+224]);
                              _mm512_store_ps(&zim[i+224], _mm512_sub_ps(zmm29,zmmx));
                              register const __m512 zmm30 = _mm512_load_ps(&xre[i+240]);
                              _mm512_store_ps(&zre[i+240], _mm512_sub_ps(zmm30,zmmx));
                              register const __m512 zmm31 = _mm512_load_ps(&xim[i+240]);
                              _mm512_storeu_ps(&zim[i+240], _mm512_sub_ps(zmm31,zmmx)); 
#else
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+255],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+255],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              register const __m512 zmm12 = _mm512_load_ps(&xre[i+96]);
                              register const __m512 zmm13 = _mm512_load_ps(&xim[i+96]);
                              register const __m512 zmm14 = _mm512_load_ps(&xre[i+112]);
                              register const __m512 zmm15 = _mm512_load_ps(&xim[i+112]);
                              register const __m512 zmm16 = _mm512_load_ps(&xre[i+128]);
                              register const __m512 zmm17 = _mm512_load_ps(&xim[i+128]);
                              register const __m512 zmm18 = _mm512_load_ps(&xre[i+144]);
                              register const __m512 zmm19 = _mm512_load_ps(&xim[i+144]);
                              register const __m512 zmm20 = _mm512_load_ps(&xre[i+160]);
                              register const __m512 zmm21 = _mm512_load_ps(&xim[i+160]);
                              register const __m512 zmm22 = _mm512_load_ps(&xre[i+176]);
                              register const __m512 zmm23 = _mm512_load_ps(&xim[i+176]);
                              register const __m512 zmm24 = _mm512_load_ps(&xre[i+192]);
                              register const __m512 zmm25 = _mm512_load_ps(&xim[i+192]);
                              register const __m512 zmm26 = _mm512_load_ps(&xre[i+208]);
                              register const __m512 zmm27 = _mm512_load_ps(&xim[i+208]);
                              register const __m512 zmm28 = _mm512_load_ps(&xre[i+224]);
                              register const __m512 zmm29 = _mm512_load_ps(&xim[i+224]);
                              register const __m512 zmm30 = _mm512_load_ps(&xre[i+240]);
                              register const __m512 zmm31 = _mm512_load_ps(&xim[i+240]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
                              _mm512_store_ps(&zre[i+160], _mm512_sub_ps(zmm20,zmmx));
                              _mm512_store_ps(&zim[i+160], _mm512_sub_ps(zmm21,zmmx));
                              _mm512_store_ps(&zre[i+176], _mm512_sub_ps(zmm22,zmmx));
                              _mm512_store_ps(&zim[i+176], _mm512_sub_ps(zmm23,zmmx));
                              _mm512_store_ps(&zre[i+192], _mm512_sub_ps(zmm24,zmmx));
                              _mm512_store_ps(&zim[i+192], _mm512_sub_ps(zmm25,zmmx));
                              _mm512_store_ps(&zre[i+208], _mm512_sub_ps(zmm26,zmmx));
                              _mm512_store_ps(&zim[i+208], _mm512_sub_ps(zmm27,zmmx));
                              _mm512_store_ps(&zre[i+224], _mm512_sub_ps(zmm28,zmmx));
                              _mm512_store_ps(&zim[i+224], _mm512_sub_ps(zmm29,zmmx));
                              _mm512_store_ps(&zre[i+240], _mm512_sub_ps(zmm30,zmmx));
                              _mm512_store_ps(&zim[i+240], _mm512_sub_ps(zmm31,zmmx)); 
                     
#endif
                         }

                        for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              register const __m512 zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              register const __m512 zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              register const __m512 zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              register const __m512 zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              register const __m512 zmm12 = _mm512_load_ps(&xre[i+96]);
                              register const __m512 zmm13 = _mm512_load_ps(&xim[i+96]);
                              register const __m512 zmm14 = _mm512_load_ps(&xre[i+112]);
                              register const __m512 zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
#endif
                        }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#endif
                       }

                      for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#endif
                     }
                    
                    
                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] - s;
                          zim[i] = zim[i] - s;
                    }
                 }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                         register const __m512 zmmx    = vs;
                         const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              register const __m512 zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              register const __m512 zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              register const __m512 zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              register const __m512 zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              register const __m512 zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              register const __m512 zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              register const __m512 zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              register const __m512 zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
#else
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              register const __m512 zmm12 = _mm512_loadu_ps(&xre[i+96]);
                              register const __m512 zmm13 = _mm512_loadu_ps(&xim[i+96]);
                              register const __m512 zmm14 = _mm512_loadu_ps(&xre[i+112]);
                              register const __m512 zmm15 = _mm512_loadu_ps(&xim[i+112]);
                              register const __m512 zmm16 = _mm512_loadu_ps(&xre[i+128]);
                              register const __m512 zmm17 = _mm512_loadu_ps(&xim[i+128]);
                              register const __m512 zmm18 = _mm512_loadu_ps(&xre[i+144]);
                              register const __m512 zmm19 = _mm512_loadu_ps(&xim[i+144]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              _mm512_storeu_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm512_storeu_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              _mm512_storeu_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              _mm512_storeu_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              _mm512_storeu_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              _mm512_storeu_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              _mm512_storeu_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
#endif
                        }

                       for(; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx)); 
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
#endif
                       }

                      for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#endif
                     }

                    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#endif
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] - s;
                          zim[i] = zim[i] - s;
                    }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                         register const __m512 zmmx    = vs;
                         const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              register const __m512 zmm12 = _mm512_load_ps(&xre[i+96]);
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              register const __m512 zmm13 = _mm512_load_ps(&xim[i+96]);
                              _mm512_storeu_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              register const __m512 zmm14 = _mm512_load_ps(&xre[i+112]);
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              register const __m512 zmm15 = _mm512_load_ps(&xim[i+112]);
                              _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              register const __m512 zmm16 = _mm512_load_ps(&xre[i+128]);
                              _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              register const __m512 zmm17 = _mm512_load_ps(&xim[i+128]);
                              _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              register const __m512 zmm18 = _mm512_load_ps(&xre[i+144]);
                              _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              register const __m512 zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
#else
                              _mm_prefetch((const char *)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+128],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              register const __m512 zmm12 = _mm512_load_ps(&xre[i+96]);
                              register const __m512 zmm13 = _mm512_load_ps(&xim[i+96]);
                              register const __m512 zmm14 = _mm512_load_ps(&xre[i+112]);
                              register const __m512 zmm15 = _mm512_load_ps(&xim[i+112]);
                              register const __m512 zmm16 = _mm512_load_ps(&xre[i+128]);
                              register const __m512 zmm17 = _mm512_load_ps(&xim[i+128]);
                              register const __m512 zmm18 = _mm512_load_ps(&xre[i+144]);
                              register const __m512 zmm19 = _mm512_load_ps(&xim[i+144]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                              _mm512_store_ps(&zre[i+96], _mm512_sub_ps(zmm12,zmmx));
                              _mm512_store_ps(&zim[i+96], _mm512_sub_ps(zmm13,zmmx));
                              _mm512_store_ps(&zre[i+112], _mm512_sub_ps(zmm14,zmmx));
                              _mm512_store_ps(&zim[i+112], _mm512_sub_ps(zmm15,zmmx));
                              _mm512_store_ps(&zre[i+128], _mm512_sub_ps(zmm16,zmmx));
                              _mm512_store_ps(&zim[i+128], _mm512_sub_ps(zmm17,zmmx));
                              _mm512_store_ps(&zre[i+144], _mm512_sub_ps(zmm18,zmmx));
                              _mm512_store_ps(&zim[i+144], _mm512_sub_ps(zmm19,zmmx));
#endif
                        }

                       for(; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx)); 
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
#endif
                       }

                      for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#endif
                     }

                    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#endif
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] - s;
                          zim[i] = zim[i] - s;
                    }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                         register const __m512 zmmx    = vs;
                         const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                            
#else
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_loadu_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_loadu_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_loadu_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_loadu_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_storeu_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_storeu_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_storeu_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
#endif
                        }

                      

                      for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_loadu_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_loadu_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_loadu_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_loadu_ps(&xim[i+48]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_storeu_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_storeu_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_storeu_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_storeu_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_loadu_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_loadu_ps(&xim[i+16]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_storeu_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_storeu_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#endif
                     }

                    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#else
                              register const __m512 zmm0 = _mm512_loadu_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_loadu_ps(&xim[i+0]);
                              _mm512_storeu_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_storeu_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#endif
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] - s;
                          zim[i] = zim[i] - s;
                    }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csub_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                         register const __m512 zmmx    = vs;
                         const float * __restrict pzmm =  (const float*)&vs[0];
                         for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_storeu_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
                            
#else
                              _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              register const __m512 zmm8 = _mm512_load_ps(&xre[i+64]);
                              register const __m512 zmm9 = _mm512_load_ps(&xim[i+64]);
                              register const __m512 zmm10 = _mm512_load_ps(&xre[i+80]);
                              register const __m512 zmm11 = _mm512_load_ps(&xim[i+80]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
                              _mm512_store_ps(&zre[i+64], _mm512_sub_ps(zmm8,zmmx));
                              _mm512_store_ps(&zim[i+64], _mm512_sub_ps(zmm9,zmmx));
                              _mm512_store_ps(&zre[i+80], _mm512_sub_ps(zmm10,zmmx));
                              _mm512_store_ps(&zim[i+80], _mm512_sub_ps(zmm11,zmmx));
#endif
                        }

                      

                      for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              register const __m512 zmm4 = _mm512_load_ps(&xre[i+32]);
                              register const __m512 zmm5 = _mm512_load_ps(&xim[i+32]);
                              register const __m512 zmm6 = _mm512_load_ps(&xre[i+48]);
                              register const __m512 zmm7 = _mm512_load_ps(&xim[i+48]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
                              _mm512_store_ps(&zre[i+32], _mm512_sub_ps(zmm4,zmmx));
                              _mm512_store_ps(&zim[i+32], _mm512_sub_ps(zmm5,zmmx));
                              _mm512_store_ps(&zre[i+48], _mm512_sub_ps(zmm6,zmmx));
                              _mm512_store_ps(&zim[i+48], _mm512_sub_ps(zmm7,zmmx));
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              register const __m512 zmm2 = _mm512_load_ps(&xre[i+16]);
                              register const __m512 zmm3 = _mm512_load_ps(&xim[i+16]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
                              _mm512_store_ps(&zre[i+16], _mm512_sub_ps(zmm2,zmmx));
                              _mm512_store_ps(&zim[i+16], _mm512_sub_ps(zmm3,zmmx));
#endif
                     }

                    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#else
                              register const __m512 zmm0 = _mm512_load_ps(&xre[i+0]);
                              register const __m512 zmm1 = _mm512_load_ps(&xim[i+0]);
                              _mm512_store_ps(&zre[i+0], _mm512_sub_ps(zmm0,zmmx));
                              _mm512_store_ps(&zim[i+0], _mm512_sub_ps(zmm1,zmmx));
#endif
                    }

                    for(; (i+0) < n; i += 1) {
                          const float s = pzmm[i];
                          zre[i] = xre[i] - s;
                          zim[i] = zim[i] - s;
                    }

                 }





  







 
 




  

      } // math

} // gms





#endif /*__GMS_CSUB_VEC_ZMM16R4_HPP__*/
