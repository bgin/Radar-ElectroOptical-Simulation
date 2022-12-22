
#ifndef __GMS_CEQ_VEC_ZMM16R4_HPP__
#define __GMS_CEQ_VEC_ZMM16R4_HPP__ 221220221003


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

    const unsigned int GMS_CEQ_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CEQ_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CEQ_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CEQ_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CEQ_VEC_ZMM16R4_MAJOR+
      100U*GMS_CEQ_VEC_ZMM16R4_MINOR+
      10U*GMS_CEQ_VEC_ZMM16R4_MICRO;
    const char * const GMS_CEQ_VEC_ZMM16R4_CREATION_DATE = "22-12-2022 10:03 AM +00200 (THR 22 DEC 2022 GMT+2)";
    const char * const GMS_CEQ_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CEQ_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CEQ_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex vectors equality comparison."

}


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"



namespace  gms {

         namespace math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void ceqv_zmm16r4_unroll_16x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  __mmask16 * __restrict eqr,
                                                  __mmask16 * __restrict eqi,
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
                        int32_t i,j;
                        j = 0;
                        for(i = 0; (i+255) < n; i += 256) {
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+16]);
                            zmm5 = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&eqr[j+1],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+16]);
                            zmm7 = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&eqi[j+1],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&eqr[j+2],
                                             _mm512_cmp_ps_mask(zmm8,zmm9,_CMP_EQ_OQ));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&eqi[j+2],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&eqr[j+3],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&eqi[j+3],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&eqr[j+4],
                                             _mm512_cmp_ps_mask(zmm16,zmm17,_CMP_EQ_OQ));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&eqi[j+4],
                                             _mm512_cmp_ps_mask(zmm18,zmm19,_CMP_EQ_OQ));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&eqr[j+5],
                                             _mm512_cmp_ps_mask(zmm20,zmm21,_CMP_EQ_OQ));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&eqi[j+5],
                                             _mm512_cmp_ps_mask(zmm22,zmm23,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&eqr[j+6],
                                             _mm512_cmp_ps_mask(zmm24,zmm25,_CMP_EQ_OQ));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&eqi[j+6],
                                             _mm512_cmp_ps_mask(zmm26,zmm27,_CMP_EQ_OQ));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&eqr[j+7],
                                             _mm512_cmp_ps_mask(zmm28,zmm29,_CMP_EQ_OQ));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&eqi[j+7],
                                             _mm512_cmp_ps_mask(zmm30,zmm31,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                            zmm0 = _mm512_loadu_ps(&xre[i+128]);
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&eqr[j+8],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&eqi[j+8],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+144]);
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&eqr[j+9],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&eqi[j+9],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+192],_MM_HINT_T0);
                            zmm8 = _mm512_loadu_ps(&xre[i+160]);
                            zmm9 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&eqr[j+10],
                                             _mm512_cmp_ps_mask(zmm9,zmm8,_CMP_EQ_OQ));
                            zmm10= _mm512_loadu_ps(&xim[i+160]);
                            zmm11= _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&eqi[j+10],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+176]);
                            zmm13 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&eqr[j+11],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+176]);
                            zmm15 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&eqi[j+11],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+224],_MM_HINT_T0);
                            zmm16 = _mm512_loadu_ps(&xre[i+192]);
                            zmm17 = _mm512_loadu_ps(&yre[i+192]);
                            _mm512_storeu_ps(&eqr[j+12],
                                             _mm512_cmp_ps_mask(zmm16,zmm17,_CMP_EQ_OQ));
                            zmm18 = _mm512_loadu_ps(&xim[i+192]);
                            zmm19 = _mm512_loadu_ps(&yim[i+192]);
                            _mm512_storeu_ps(&eqi[j+12],
                                             _mm512_cmp_ps_mask(zmm18,zmm19,_CMP_EQ_OQ));
                            zmm20 = _mm512_loadu_ps(&xre[i+208]);
                            zmm21 = _mm512_loadu_ps(&yre[i+208]);
                            _mm512_storeu_ps(&eqr[j+13],
                                             _mm512_cmp_ps_mask(zmm20,zmm21,_CMP_EQ_OQ));
                            zmm22 = _mm512_loadu_ps(&xim[i+208]);
                            zmm23 = _mm512_loadu_ps(&yim[i+208]);
                            _mm512_storeu_ps(&eqi[j+13],
                                             _mm512_cmp_ps_mask(zmm22,zmm23,_CMP_EQ_OQ));
                            _mm_prefetch((const char*)&xre[i+256],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yre[i+256],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+256],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+256],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&xre[i+224]);
                            zmm25 = _mm512_loadu_ps(&yre[i+224]);
                            _mm512_storeu_ps(&eqr[j+14],
                                             _mm512_cmp_ps_mask(zmm24,zmm25,_CMP_EQ_OQ));
                            zmm26 = _mm512_loadu_ps(&xim[i+224]);
                            zmm27 = _mm512_loadu_ps(&yim[i+224]);
                            _mm512_storeu_ps(&eqi[j+14],
                                             _mm512_cmp_ps_mask(zmm26,zmm27,_CMP_EQ_OQ));
                            zmm28 = _mm512_loadu_ps(&xre[i+240]);
                            zmm29 = _mm512_loadu_ps(&yre[i+240]);
                            _mm512_storeu_ps(&eqr[j+15],
                                             _mm512_cmp_ps_mask(zmm28,zmm29,_CMP_EQ_OQ));
                            zmm30 = _mm512_loadu_ps(&xim[i+240]);
                            zmm31 = _mm512_loadu_ps(&yim[i+240]);
                            _mm512_storeu_ps(&eqi[j+15],
                                             _mm512_cmp_ps_mask(zmm30,zmm31,_CMP_EQ_OQ));
                            j += 16;
                       }

                        for(; (i+191) < n; i += 192) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+16]);
                            zmm5 = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&eqr[j+1],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+16]);
                            zmm7 = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&eqi[j+1],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&eqr[j+2],
                                             _mm512_cmp_ps_mask(zmm8,zmm9,_CMP_EQ_OQ));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&eqi[j+2],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&eqr[j+3],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&eqi[j+3],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&eqr[j+4],
                                             _mm512_cmp_ps_mask(zmm16,zmm17,_CMP_EQ_OQ));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&eqi[j+4],
                                             _mm512_cmp_ps_mask(zmm18,zmm19,_CMP_EQ_OQ));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&eqr[j+5],
                                             _mm512_cmp_ps_mask(zmm20,zmm21,_CMP_EQ_OQ));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&eqi[j+5],
                                             _mm512_cmp_ps_mask(zmm22,zmm23,_CMP_EQ_OQ));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&eqr[j+6],
                                             _mm512_cmp_ps_mask(zmm24,zmm25,_CMP_EQ_OQ));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&eqi[j+6],
                                             _mm512_cmp_ps_mask(zmm26,zmm27,_CMP_EQ_OQ));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&eqr[j+7],
                                             _mm512_cmp_ps_mask(zmm28,zmm29,_CMP_EQ_OQ));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&eqi[j+7],
                                             _mm512_cmp_ps_mask(zmm30,zmm31,_CMP_EQ_OQ));
                            zmm0 = _mm512_loadu_ps(&xre[i+128]);
                            zmm1 = _mm512_loadu_ps(&yre[i+128]);
                            _mm512_storeu_ps(&eqr[j+8],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+128]);
                            zmm3 = _mm512_loadu_ps(&yim[i+128]);
                            _mm512_storeu_ps(&eqi[j+8],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+144]);
                            zmm5 = _mm512_loadu_ps(&yre[i+144]);
                            _mm512_storeu_ps(&eqr[j+9],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+144]);
                            zmm7 = _mm512_loadu_ps(&yim[i+144]);
                            _mm512_storeu_ps(&eqi[j+9],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            zmm8 = _mm512_loadu_ps(&xre[i+160]);
                            zmm9 = _mm512_loadu_ps(&yre[i+160]);
                            _mm512_storeu_ps(&eqr[j+10],
                                             _mm512_cmp_ps_mask(zmm9,zmm8,_CMP_EQ_OQ));
                            zmm10= _mm512_loadu_ps(&xim[i+160]);
                            zmm11= _mm512_loadu_ps(&yim[i+160]);
                            _mm512_storeu_ps(&eqi[j+10],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+176]);
                            zmm13 = _mm512_loadu_ps(&yre[i+176]);
                            _mm512_storeu_ps(&eqr[j+11],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+176]);
                            zmm15 = _mm512_loadu_ps(&yim[i+176]);
                            _mm512_storeu_ps(&eqi[j+11],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            j += 16;
                       }

                        for(; (i+127) < n; i += 128) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+16]);
                            zmm5 = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&eqr[j+1],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+16]);
                            zmm7 = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&eqi[j+1],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&eqr[j+2],
                                             _mm512_cmp_ps_mask(zmm8,zmm9,_CMP_EQ_OQ));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&eqi[j+2],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&eqr[j+3],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&eqi[j+3],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&eqr[j+4],
                                             _mm512_cmp_ps_mask(zmm16,zmm17,_CMP_EQ_OQ));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&eqi[j+4],
                                             _mm512_cmp_ps_mask(zmm18,zmm19,_CMP_EQ_OQ));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&yre[i+80]);
                            _mm512_storeu_ps(&eqr[j+5],
                                             _mm512_cmp_ps_mask(zmm20,zmm21,_CMP_EQ_OQ));
                            zmm22 = _mm512_loadu_ps(&xim[i+80]);
                            zmm23 = _mm512_loadu_ps(&yim[i+80]);
                            _mm512_storeu_ps(&eqi[j+5],
                                             _mm512_cmp_ps_mask(zmm22,zmm23,_CMP_EQ_OQ));
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&yre[i+96]);
                            _mm512_storeu_ps(&eqr[j+6],
                                             _mm512_cmp_ps_mask(zmm24,zmm25,_CMP_EQ_OQ));
                            zmm26 = _mm512_loadu_ps(&xim[i+96]);
                            zmm27 = _mm512_loadu_ps(&yim[i+96]);
                            _mm512_storeu_ps(&eqi[j+6],
                                             _mm512_cmp_ps_mask(zmm26,zmm27,_CMP_EQ_OQ));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&yre[i+112]);
                            _mm512_storeu_ps(&eqr[j+7],
                                             _mm512_cmp_ps_mask(zmm28,zmm29,_CMP_EQ_OQ));
                            zmm30 = _mm512_loadu_ps(&xim[i+112]);
                            zmm31 = _mm512_loadu_ps(&yim[i+112]);
                            _mm512_storeu_ps(&eqi[j+7],
                                             _mm512_cmp_ps_mask(zmm30,zmm31,_CMP_EQ_OQ));
                            j += 16;
                       }

                        for(; (i+79) < n; i += 80) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+16]);
                            zmm5 = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&eqr[j+1],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+16]);
                            zmm7 = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&eqi[j+1],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&eqr[j+2],
                                             _mm512_cmp_ps_mask(zmm8,zmm9,_CMP_EQ_OQ));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&eqi[j+2],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&eqr[j+3],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&eqi[j+3],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&yre[i+64]);
                            _mm512_storeu_ps(&eqr[j+4],
                                             _mm512_cmp_ps_mask(zmm16,zmm17,_CMP_EQ_OQ));
                            zmm18 = _mm512_loadu_ps(&xim[i+64]);
                            zmm19 = _mm512_loadu_ps(&yim[i+64]);
                            _mm512_storeu_ps(&eqi[j+4],
                                             _mm512_cmp_ps_mask(zmm18,zmm19,_CMP_EQ_OQ));
                            j += 16;
                       }

                         for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+16]);
                            zmm5 = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&eqr[j+1],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+16]);
                            zmm7 = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&eqi[j+1],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            zmm8 = _mm512_loadu_ps(&xre[i+32]);
                            zmm9 = _mm512_loadu_ps(&yre[i+32]);
                            _mm512_storeu_ps(&eqr[j+2],
                                             _mm512_cmp_ps_mask(zmm8,zmm9,_CMP_EQ_OQ));
                            zmm10 = _mm512_loadu_ps(&xim[i+32]);
                            zmm11 = _mm512_loadu_ps(&yim[i+32]);
                            _mm512_storeu_ps(&eqi[j+2],
                                             _mm512_cmp_ps_mask(zmm10,zmm11,_CMP_EQ_OQ));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&yre[i+48]);
                            _mm512_storeu_ps(&eqr[j+3],
                                             _mm512_cmp_ps_mask(zmm12,zmm13,_CMP_EQ_OQ));
                            zmm14 = _mm512_loadu_ps(&xim[i+48]);
                            zmm15 = _mm512_loadu_ps(&yim[i+48]);
                            _mm512_storeu_ps(&eqi[j+3],
                                             _mm512_cmp_ps_mask(zmm14,zmm15,_CMP_EQ_OQ));
                            j += 16;
                       }

                         for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            zmm4 = _mm512_loadu_ps(&xre[i+16]);
                            zmm5 = _mm512_loadu_ps(&yre[i+16]);
                            _mm512_storeu_ps(&eqr[j+1],
                                             _mm512_cmp_ps_mask(zmm4,zmm5,_CMP_EQ_OQ));
                            zmm6 = _mm512_loadu_ps(&xim[i+16]);
                            zmm7 = _mm512_loadu_ps(&yim[i+16]);
                            _mm512_storeu_ps(&eqi[j+1],
                                             _mm512_cmp_ps_mask(zmm6,zmm7,_CMP_EQ_OQ));
                            j += 16;
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&yre[i+0]);
                            _mm512_storeu_ps(&eqr[j+0],
                                             _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                            zmm2 = _mm512_loadu_ps(&xim[i+0]);
                            zmm3 = _mm512_loadu_ps(&yim[i+0]);
                            _mm512_storeu_ps(&eqi[j+0],
                                             _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
                            j += 15;
                       }
                        return;
                      // The length of 'n' shall not cause slip down to scalar loop.
                 }


                 



      } // math

} // gms












#endif /*__GMS_CEQ_VEC_ZMM16R4_HPP__*/
