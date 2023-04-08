

#ifndef __GMS_CDIV_SMITH_VEC_ZMM16R4_HPP__
#define __GMS_CDIV_SMITH_VEC_ZMM16R4_HPP__


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

    const unsigned int GMS_CDIV_SMITH_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CDIV_SMITH_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CDIV_SMITH_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CDIV_SMITH_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CDIV_SMITH_VEC_ZMM16R4_MAJOR+
      100U*GMS_CDIV_SMITH_VEC_ZMM16R4_MINOR+
      10U*GMS_CDIV_SMITH_VEC_ZMM16R4_MICRO;
    const char * const GMS_CDIV_SMITH_VEC_ZMM16R4_CREATION_DATE = "08-04-2023 09:51 AM +00200 ( SAT 08 APR 2023 GMT+2)";
    const char * const GMS_CDIV_SMITH_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CDIV_SMITH_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CDIV_SMITH_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector smith division operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_cephes.h"

namespace  gms {


        namespace math {

                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void  cdivv_smith_zmm16r4_u10x_u( const float * __restrict  xre,
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
                             _mm_prefetch((const char *)&xre[i+48],_MM_HINT_T0);
                             _mm_prefetch((const char *)&yre[i+48],_MM_HINT_T0);
                             _mm_prefetch((const char *)&xim[i+48],_MM_HINT_T0);
                             _mm_prefetch((const char *)&yim[i+48],_MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&yre[i+0]); // c
                             zmm1 = _mm512_loadu_ps(&yim[i+0]); // d
                             zmm2 = _mm512_loadu_ps(&xre[i+0]); // a
                             zmm3 = _mm512_loadu_ps(&xim[i+0]); // b
                             const __mmask16 m0 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                            _mm512_abs_ps(zmm1),
                                                            _CMP_GE_OQ);
                             zmm4 = _mm512_mask_blend_ps(m0,_mm512_div_ps(zmm0,zmm1),
                                                            _mm512_div_ps(zmm1,zmm0)); // r = zmm4
                             zmm5 = _mm512_mask_blend_ps(m0,_mm512_fmadd_ps(zmm4,zmm0,zmm1),
                                                            _mm512_fmadd_ps(zmm4,zmm1,zmm0)); // den = zmm5
                             _mm512_storeu_ps(&zre[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm2,zmm4,zmm3),zmm5),
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm3,zmm4,zmm2),zmm5)));
                             _mm512_storeu_ps(&zim[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmsub_ps(zmm3,zmm4,zmm2),zmm5),
                                                            _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,zmm4)),zmm5)));
                             zmm6 = _mm512_loadu_ps(&yre[i+16]); // c
                             zmm7 = _mm512_loadu_ps(&yim[i+16]); // d
                             zmm8 = _mm512_loadu_ps(&xre[i+16]); // a
                             zmm9 = _mm512_loadu_ps(&xim[i+16]); // b
                             const __mmask16 m1 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm6),
                                                              _mm512_abs_ps(zmm7),
                                                              _CMP_GE_OQ);
                             zmm10 = _mm512_mask_blend_ps(m1,_mm512_div_ps(zmm6,zmm7),
                                                             _mm512_div_ps(zmm7,zmm6)); // r = zmm10
                             zmm11 = _mm512_mask_blend_ps(m1,_mm512_fmadd_ps(zmm10,zmm6,zmm7),
                                                             _mm512_fmadd_ps(zmm10,zmm7,zmm6)); // den = zmm11
                             _mm512_storeu_ps(&zre[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm8,zmm10,zmm9),zmm11),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm9,zmm10,zmm8),zmm11)));
                             _mm512_storeu_ps(&zim[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm9,zmm10,zmm8),zmm11),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm9,_mm512_mul_ps(zmm8,zmm10)),zmm11)));
                             zmm12 = _mm512_loadu_ps(&yre[i+32]); // c
                             zmm13 = _mm512_loadu_ps(&yim[i+32]); // d
                             zmm14 = _mm512_loadu_ps(&xre[i+32]); // a
                             zmm15 = _mm512_loadu_ps(&xim[i+32]); // b
                             const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm12),
                                                              _mm512_abs_ps(zmm13),
                                                              _CMP_GE_OQ);
                             zmm16 = _mm512_mask_blend_ps(m2,_mm512_div_ps(zmm12,zmm13),
                                                             _mm512_div_ps(zmm13,zmm12)); // r = zmm16
                             zmm17 = _mm512_mask_blend_ps(m2,_mm512_fmadd_ps(zmm16,zmm12,zmm13),
                                                             _mm512_fmadd_ps(zmm16,zmm13,zmm12)); // den = zmm17
                             _mm512_storeu_ps(&zre[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm14,zmm16,zmm15),zmm17),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm15,zmm16,zmm14),zmm17)));
                             _mm512_storeu_ps(&zim[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm15,zmm16,zmm14),zmm17),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm15,_mm512_mul_ps(zmm14,zmm16)),zmm17)));
                             _mm_prefetch((const char *)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char *)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char *)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char *)&yim[i+96],_MM_HINT_T0);
                             zmm18 = _mm512_loadu_ps(&yre[i+48]); // c
                             zmm19 = _mm512_loadu_ps(&yim[i+48]); // d
                             zmm20 = _mm512_loadu_ps(&xre[i+48]); // a
                             zmm21 = _mm512_loadu_ps(&xim[i+48]); // b
                             const __mmask16 m3 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm18),
                                                              _mm512_abs_ps(zmm19),
                                                              _CMP_GE_OQ);
                             zmm22 = _mm512_mask_blend_ps(m3,_mm512_div_ps(zmm18,zmm19),
                                                             _mm512_div_ps(zmm19,zmm18)); // r = zmm22
                             zmm23 = _mm512_mask_blend_ps(m3,_mm512_fmadd_ps(zmm22,zmm18,zmm19),
                                                             _mm512_fmadd_ps(zmm22,zmm19,zmm18)); // den = zmm23
                             _mm512_storeu_ps(&zre[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm20,zmm22,zmm21),zmm23),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm21,zmm22,zmm20),zmm23)));
                             _mm512_storeu_ps(&zim[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm21,zmm22,zmm20),zmm23),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm21,_mm512_mul_ps(zmm20,zmm22)),zmm23)));
                             zmm24 = _mm512_loadu_ps(&yre[i+64]); // c
                             zmm25 = _mm512_loadu_ps(&yim[i+64]); // d
                             zmm26 = _mm512_loadu_ps(&xre[i+64]); // a
                             zmm27 = _mm512_loadu_ps(&xim[i+64]); // b
                             const __mmask16 m4 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm24),
                                                              _mm512_abs_ps(zmm25),
                                                              _CMP_GE_OQ);
                             zmm28 = _mm512_mask_blend_ps(m4,_mm512_div_ps(zmm24,zmm25),
                                                             _mm512_div_ps(zmm25,zmm24)); // r = zmm28
                             zmm29 = _mm512_mask_blend_ps(m4,_mm512_fmadd_ps(zmm28,zmm24,zmm25),
                                                             _mm512_fmadd_ps(zmm28,zmm25,zmm24)); // den = zmm29
                             _mm512_storeu_ps(&zre[i+64], _mm512_mask_blend_ps(m4,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm26,zmm28,zmm27),zmm29),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm27,zmm28,zmm26),zmm29)));
                             _mm512_storeu_ps(&zim[i+64], _mm512_mask_blend_ps(m4,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm27,zmm28,zmm26),zmm29),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm27,_mm512_mul_ps(zmm26,zmm28)),zmm29)));
                             zmm30 = _mm512_loadu_ps(&yre[i+80]); // c
                             zmm31 = _mm512_loadu_ps(&yim[i+80]); // d
                             zmm0  = _mm512_loadu_ps(&xre[i+80]); // a
                             zmm1  = _mm512_loadu_ps(&xim[i+80]); // b
                             const __mmask16 m5 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm30),
                                                              _mm512_abs_ps(zmm31),
                                                              _CMP_GE_OQ);
                             zmm2  = _mm512_mask_blend_ps(m5,_mm512_div_ps(zmm30,zmm31),
                                                             _mm512_div_ps(zmm31,zmm30)); // r = zmm2
                             zmm3  = _mm512_mask_blend_ps(m5,_mm512_fmadd_ps(zmm2,zmm30,zmm31),
                                                             _mm512_fmadd_ps(zmm2,zmm31,zmm30)); // den = zmm3
                             _mm512_storeu_ps(&zre[i+80], _mm512_mask_blend_ps(m5,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm0,zmm2,zmm1),zmm3),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm1,zmm2,zmm0),zmm3)));
                             _mm512_storeu_ps(&zim[i+80], _mm512_mask_blend_ps(m5,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm1,zmm2,zmm0),zmm3),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm1,_mm512_mul_ps(zmm0,zmm2)),zmm3)));
                             _mm_prefetch((const char *)&xre[i+144],_MM_HINT_T0);
                             _mm_prefetch((const char *)&yre[i+144],_MM_HINT_T0);
                             _mm_prefetch((const char *)&xim[i+144],_MM_HINT_T0);
                             _mm_prefetch((const char *)&yim[i+144],_MM_HINT_T0);
                             zmm4 = _mm512_loadu_ps(&yre[i+96]); // c
                             zmm5 = _mm512_loadu_ps(&yim[i+96]); // d
                             zmm6 = _mm512_loadu_ps(&xre[i+96]); // a
                             zmm7 = _mm512_loadu_ps(&xim[i+96]); // b
                             const __mmask16 m6 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm4),
                                                              _mm512_abs_ps(zmm5),
                                                              _CMP_GE_OQ);
                             zmm8 = _mm512_mask_blend_ps(m6,_mm512_div_ps(zmm4,zmm5),
                                                             _mm512_div_ps(zmm5,zmm4)); // r = zmm8
                             zmm9 = _mm512_mask_blend_ps(m6,_mm512_fmadd_ps(zmm8,zmm4,zmm5),
                                                             _mm512_fmadd_ps(zmm8,zmm5,zmm4)); // den = zmm9
                             _mm512_storeu_ps(&zre[i+96], _mm512_mask_blend_ps(m6,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm6,zmm8,zmm7),zmm9),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm7,zmm8,zmm6),zmm9)));
                             _mm512_storeu_ps(&zim[i+96], _mm512_mask_blend_ps(m6,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm7,zmm8,zmm6),zmm9),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm7,_mm512_mul_ps(zmm6,zmm8)),zmm9)));
                             zmm10 = _mm512_loadu_ps(&yre[i+112]); // c
                             zmm11 = _mm512_loadu_ps(&yim[i+112]); // d
                             zmm12 = _mm512_loadu_ps(&xre[i+112]); // a
                             zmm13 = _mm512_loadu_ps(&xim[i+112]); // b
                             const __mmask16 m7 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm10),
                                                              _mm512_abs_ps(zmm11),
                                                              _CMP_GE_OQ);
                             zmm14 = _mm512_mask_blend_ps(m7,_mm512_div_ps(zmm10,zmm11),
                                                             _mm512_div_ps(zmm11,zmm10)); // r = zmm14
                             zmm15 = _mm512_mask_blend_ps(m7,_mm512_fmadd_ps(zmm14,zmm10,zmm11),
                                                             _mm512_fmadd_ps(zmm14,zmm11,zmm10)); // den = zmm15
                             _mm512_storeu_ps(&zre[i+112], _mm512_mask_blend_ps(m7,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm12,zmm14,zmm13),zmm15),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm13,zmm14,zmm12),zmm15)));
                             _mm512_storeu_ps(&zim[i+112], _mm512_mask_blend_ps(m7,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm13,zmm14,zmm12),zmm15),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm13,_mm512_mul_ps(zmm12,zmm14)),zmm15)));
                             zmm16 = _mm512_loadu_ps(&yre[i+128]); // c
                             zmm17 = _mm512_loadu_ps(&yim[i+128]); // d
                             zmm18 = _mm512_loadu_ps(&xre[i+128]); // a
                             zmm19 = _mm512_loadu_ps(&xim[i+128]); // b
                             const __mmask16 m8 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm16),
                                                              _mm512_abs_ps(zmm17),
                                                              _CMP_GE_OQ);
                             zmm20 = _mm512_mask_blend_ps(m8,_mm512_div_ps(zmm16,zmm17),
                                                             _mm512_div_ps(zmm17,zmm16)); // r = zmm20
                             zmm21 = _mm512_mask_blend_ps(m8,_mm512_fmadd_ps(zmm20,zmm16,zmm17),
                                                             _mm512_fmadd_ps(zmm20,zmm17,zmm16)); // den = zmm21
                             _mm512_storeu_ps(&zre[i+128], _mm512_mask_blend_ps(m8,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm18,zmm20,zmm19),zmm21),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm19,zmm20,zmm18),zmm21)));
                             _mm512_storeu_ps(&zim[i+128], _mm512_mask_blend_ps(m8,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm19,zmm20,zmm18),zmm21),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm19,_mm512_mul_ps(zmm18,zmm20)),zmm21)));
                             zmm22 = _mm512_loadu_ps(&yre[i+144]); // c
                             zmm23 = _mm512_loadu_ps(&yim[i+144]); // d
                             zmm24 = _mm512_loadu_ps(&xre[i+144]); // a
                             zmm25 = _mm512_loadu_ps(&xim[i+144]); // b
                             const __mmask18 m9 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm22),
                                                              _mm512_abs_ps(zmm23),
                                                              _CMP_GE_OQ);
                             zmm26 = _mm512_mask_blend_ps(m9,_mm512_div_ps(zmm22,zmm23),
                                                             _mm512_div_ps(zmm23,zmm22)); // r = zmm26
                             zmm27 = _mm512_mask_blend_ps(m9,_mm512_fmadd_ps(zmm26,zmm22,zmm23),
                                                             _mm512_fmadd_ps(zmm26,zmm23,zmm22)); // den = zmm27
                             _mm512_storeu_ps(&zre[i+144], _mm512_mask_blend_ps(m9,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm24,zmm26,zmm25),zmm27),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm25,zmm26,zmm24),zmm27)));
                             _mm512_storeu_ps(&zim[i+144], _mm512_mask_blend_ps(m9,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm25,zmm26,zmm24),zmm27),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm25,_mm512_mul_ps(zmm24,zmm26)),zmm27)));
                        }

                          for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_loadu_ps(&yre[i+0]); // c
                             zmm1 = _mm512_loadu_ps(&yim[i+0]); // d
                             zmm2 = _mm512_loadu_ps(&xre[i+0]); // a
                             zmm3 = _mm512_loadu_ps(&xim[i+0]); // b
                             const __mmask16 m0 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                            _mm512_abs_ps(zmm1),
                                                            _CMP_GE_OQ);
                             zmm4 = _mm512_mask_blend_ps(m0,_mm512_div_ps(zmm0,zmm1),
                                                            _mm512_div_ps(zmm1,zmm0)); // r = zmm4
                             zmm5 = _mm512_mask_blend_ps(m0,_mm512_fmadd_ps(zmm4,zmm0,zmm1),
                                                            _mm512_fmadd_ps(zmm4,zmm1,zmm0)); // den = zmm5
                             _mm512_storeu_ps(&zre[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm2,zmm4,zmm3),zmm5),
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm3,zmm4,zmm2),zmm5)));
                             _mm512_storeu_ps(&zim[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmsub_ps(zmm3,zmm4,zmm2),zmm5),
                                                            _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,zmm4)),zmm5)));
                             zmm6 = _mm512_loadu_ps(&yre[i+16]); // c
                             zmm7 = _mm512_loadu_ps(&yim[i+16]); // d
                             zmm8 = _mm512_loadu_ps(&xre[i+16]); // a
                             zmm9 = _mm512_loadu_ps(&xim[i+16]); // b
                             const __mmask16 m1 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm6),
                                                              _mm512_abs_ps(zmm7),
                                                              _CMP_GE_OQ);
                             zmm10 = _mm512_mask_blend_ps(m1,_mm512_div_ps(zmm6,zmm7),
                                                             _mm512_div_ps(zmm7,zmm6)); // r = zmm10
                             zmm11 = _mm512_mask_blend_ps(m1,_mm512_fmadd_ps(zmm10,zmm6,zmm7),
                                                             _mm512_fmadd_ps(zmm10,zmm7,zmm6)); // den = zmm11
                             _mm512_storeu_ps(&zre[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm8,zmm10,zmm9),zmm11),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm9,zmm10,zmm8),zmm11)));
                             _mm512_storeu_ps(&zim[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm9,zmm10,zmm8),zmm11),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm9,_mm512_mul_ps(zmm8,zmm10)),zmm11)));
                             zmm12 = _mm512_loadu_ps(&yre[i+32]); // c
                             zmm13 = _mm512_loadu_ps(&yim[i+32]); // d
                             zmm14 = _mm512_loadu_ps(&xre[i+32]); // a
                             zmm15 = _mm512_loadu_ps(&xim[i+32]); // b
                             const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm12),
                                                              _mm512_abs_ps(zmm13),
                                                              _CMP_GE_OQ);
                             zmm16 = _mm512_mask_blend_ps(m2,_mm512_div_ps(zmm12,zmm13),
                                                             _mm512_div_ps(zmm13,zmm12)); // r = zmm16
                             zmm17 = _mm512_mask_blend_ps(m2,_mm512_fmadd_ps(zmm16,zmm12,zmm13),
                                                             _mm512_fmadd_ps(zmm16,zmm13,zmm12)); // den = zmm17
                             _mm512_storeu_ps(&zre[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm14,zmm16,zmm15),zmm17),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm15,zmm16,zmm14),zmm17)));
                             _mm512_storeu_ps(&zim[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm15,zmm16,zmm14),zmm17),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm15,_mm512_mul_ps(zmm14,zmm16)),zmm17)));
                             zmm18 = _mm512_loadu_ps(&yre[i+48]); // c
                             zmm19 = _mm512_loadu_ps(&yim[i+48]); // d
                             zmm20 = _mm512_loadu_ps(&xre[i+48]); // a
                             zmm21 = _mm512_loadu_ps(&xim[i+48]); // b
                             const __mmask16 m3 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm18),
                                                              _mm512_abs_ps(zmm19),
                                                              _CMP_GE_OQ);
                             zmm22 = _mm512_mask_blend_ps(m3,_mm512_div_ps(zmm18,zmm19),
                                                             _mm512_div_ps(zmm19,zmm18)); // r = zmm22
                             zmm23 = _mm512_mask_blend_ps(m3,_mm512_fmadd_ps(zmm22,zmm18,zmm19),
                                                             _mm512_fmadd_ps(zmm22,zmm19,zmm18)); // den = zmm23
                             _mm512_storeu_ps(&zre[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm20,zmm22,zmm21),zmm23),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm21,zmm22,zmm20),zmm23)));
                             _mm512_storeu_ps(&zim[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm21,zmm22,zmm20),zmm23),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm21,_mm512_mul_ps(zmm20,zmm22)),zmm23)));
                             zmm24 = _mm512_loadu_ps(&yre[i+64]); // c
                             zmm25 = _mm512_loadu_ps(&yim[i+64]); // d
                             zmm26 = _mm512_loadu_ps(&xre[i+64]); // a
                             zmm27 = _mm512_loadu_ps(&xim[i+64]); // b
                             const __mmask16 m4 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm24),
                                                              _mm512_abs_ps(zmm25),
                                                              _CMP_GE_OQ);
                             zmm28 = _mm512_mask_blend_ps(m4,_mm512_div_ps(zmm24,zmm25),
                                                             _mm512_div_ps(zmm25,zmm24)); // r = zmm28
                             zmm29 = _mm512_mask_blend_ps(m4,_mm512_fmadd_ps(zmm28,zmm24,zmm25),
                                                             _mm512_fmadd_ps(zmm28,zmm25,zmm24)); // den = zmm29
                             _mm512_storeu_ps(&zre[i+64], _mm512_mask_blend_ps(m4,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm26,zmm28,zmm27),zmm29),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm27,zmm28,zmm26),zmm29)));
                             _mm512_storeu_ps(&zim[i+64], _mm512_mask_blend_ps(m4,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm27,zmm28,zmm26),zmm29),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm27,_mm512_mul_ps(zmm26,zmm28)),zmm29)));
                             zmm30 = _mm512_loadu_ps(&yre[i+80]); // c
                             zmm31 = _mm512_loadu_ps(&yim[i+80]); // d
                             zmm0  = _mm512_loadu_ps(&xre[i+80]); // a
                             zmm1  = _mm512_loadu_ps(&xim[i+80]); // b
                             const __mmask16 m5 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm30),
                                                              _mm512_abs_ps(zmm31),
                                                              _CMP_GE_OQ);
                             zmm2  = _mm512_mask_blend_ps(m5,_mm512_div_ps(zmm30,zmm31),
                                                             _mm512_div_ps(zmm31,zmm30)); // r = zmm2
                             zmm3  = _mm512_mask_blend_ps(m5,_mm512_fmadd_ps(zmm2,zmm30,zmm31),
                                                             _mm512_fmadd_ps(zmm2,zmm31,zmm30)); // den = zmm3
                             _mm512_storeu_ps(&zre[i+80], _mm512_mask_blend_ps(m5,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm0,zmm2,zmm1),zmm3),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm1,zmm2,zmm0),zmm3)));
                             _mm512_storeu_ps(&zim[i+80], _mm512_mask_blend_ps(m5,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm1,zmm2,zmm0),zmm3),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm1,_mm512_mul_ps(zmm0,zmm2)),zmm3)));
                             zmm4 = _mm512_loadu_ps(&yre[i+96]); // c
                             zmm5 = _mm512_loadu_ps(&yim[i+96]); // d
                             zmm6 = _mm512_loadu_ps(&xre[i+96]); // a
                             zmm7 = _mm512_loadu_ps(&xim[i+96]); // b
                             const __mmask16 m6 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm4),
                                                              _mm512_abs_ps(zmm5),
                                                              _CMP_GE_OQ);
                             zmm8 = _mm512_mask_blend_ps(m6,_mm512_div_ps(zmm4,zmm5),
                                                             _mm512_div_ps(zmm5,zmm4)); // r = zmm8
                             zmm9 = _mm512_mask_blend_ps(m6,_mm512_fmadd_ps(zmm8,zmm4,zmm5),
                                                             _mm512_fmadd_ps(zmm8,zmm5,zmm4)); // den = zmm9
                             _mm512_storeu_ps(&zre[i+96], _mm512_mask_blend_ps(m6,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm6,zmm8,zmm7),zmm9),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm7,zmm8,zmm6),zmm9)));
                             _mm512_storeu_ps(&zim[i+96], _mm512_mask_blend_ps(m6,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm7,zmm8,zmm6),zmm9),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm7,_mm512_mul_ps(zmm6,zmm8)),zmm9)));
                             zmm10 = _mm512_loadu_ps(&yre[i+112]); // c
                             zmm11 = _mm512_loadu_ps(&yim[i+112]); // d
                             zmm12 = _mm512_loadu_ps(&xre[i+112]); // a
                             zmm13 = _mm512_loadu_ps(&xim[i+112]); // b
                             const __mmask16 m7 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm10),
                                                              _mm512_abs_ps(zmm11),
                                                              _CMP_GE_OQ);
                             zmm14 = _mm512_mask_blend_ps(m7,_mm512_div_ps(zmm10,zmm11),
                                                             _mm512_div_ps(zmm11,zmm10)); // r = zmm14
                             zmm15 = _mm512_mask_blend_ps(m7,_mm512_fmadd_ps(zmm14,zmm10,zmm11),
                                                             _mm512_fmadd_ps(zmm14,zmm11,zmm10)); // den = zmm15
                             _mm512_storeu_ps(&zre[i+112], _mm512_mask_blend_ps(m7,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm12,zmm14,zmm13),zmm15),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm13,zmm14,zmm12),zmm15)));
                             _mm512_storeu_ps(&zim[i+112], _mm512_mask_blend_ps(m7,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm13,zmm14,zmm12),zmm15),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm13,_mm512_mul_ps(zmm12,zmm14)),zmm15)));
                        }

                          for(; (i+95) < n; i += 96) {
                             zmm0 = _mm512_loadu_ps(&yre[i+0]); // c
                             zmm1 = _mm512_loadu_ps(&yim[i+0]); // d
                             zmm2 = _mm512_loadu_ps(&xre[i+0]); // a
                             zmm3 = _mm512_loadu_ps(&xim[i+0]); // b
                             const __mmask16 m0 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                            _mm512_abs_ps(zmm1),
                                                            _CMP_GE_OQ);
                             zmm4 = _mm512_mask_blend_ps(m0,_mm512_div_ps(zmm0,zmm1),
                                                            _mm512_div_ps(zmm1,zmm0)); // r = zmm4
                             zmm5 = _mm512_mask_blend_ps(m0,_mm512_fmadd_ps(zmm4,zmm0,zmm1),
                                                            _mm512_fmadd_ps(zmm4,zmm1,zmm0)); // den = zmm5
                             _mm512_storeu_ps(&zre[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm2,zmm4,zmm3),zmm5),
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm3,zmm4,zmm2),zmm5)));
                             _mm512_storeu_ps(&zim[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmsub_ps(zmm3,zmm4,zmm2),zmm5),
                                                            _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,zmm4)),zmm5)));
                             zmm6 = _mm512_loadu_ps(&yre[i+16]); // c
                             zmm7 = _mm512_loadu_ps(&yim[i+16]); // d
                             zmm8 = _mm512_loadu_ps(&xre[i+16]); // a
                             zmm9 = _mm512_loadu_ps(&xim[i+16]); // b
                             const __mmask16 m1 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm6),
                                                              _mm512_abs_ps(zmm7),
                                                              _CMP_GE_OQ);
                             zmm10 = _mm512_mask_blend_ps(m1,_mm512_div_ps(zmm6,zmm7),
                                                             _mm512_div_ps(zmm7,zmm6)); // r = zmm10
                             zmm11 = _mm512_mask_blend_ps(m1,_mm512_fmadd_ps(zmm10,zmm6,zmm7),
                                                             _mm512_fmadd_ps(zmm10,zmm7,zmm6)); // den = zmm11
                             _mm512_storeu_ps(&zre[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm8,zmm10,zmm9),zmm11),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm9,zmm10,zmm8),zmm11)));
                             _mm512_storeu_ps(&zim[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm9,zmm10,zmm8),zmm11),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm9,_mm512_mul_ps(zmm8,zmm10)),zmm11)));
                             zmm12 = _mm512_loadu_ps(&yre[i+32]); // c
                             zmm13 = _mm512_loadu_ps(&yim[i+32]); // d
                             zmm14 = _mm512_loadu_ps(&xre[i+32]); // a
                             zmm15 = _mm512_loadu_ps(&xim[i+32]); // b
                             const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm12),
                                                              _mm512_abs_ps(zmm13),
                                                              _CMP_GE_OQ);
                             zmm16 = _mm512_mask_blend_ps(m2,_mm512_div_ps(zmm12,zmm13),
                                                             _mm512_div_ps(zmm13,zmm12)); // r = zmm16
                             zmm17 = _mm512_mask_blend_ps(m2,_mm512_fmadd_ps(zmm16,zmm12,zmm13),
                                                             _mm512_fmadd_ps(zmm16,zmm13,zmm12)); // den = zmm17
                             _mm512_storeu_ps(&zre[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm14,zmm16,zmm15),zmm17),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm15,zmm16,zmm14),zmm17)));
                             _mm512_storeu_ps(&zim[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm15,zmm16,zmm14),zmm17),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm15,_mm512_mul_ps(zmm14,zmm16)),zmm17)));
                             zmm18 = _mm512_loadu_ps(&yre[i+48]); // c
                             zmm19 = _mm512_loadu_ps(&yim[i+48]); // d
                             zmm20 = _mm512_loadu_ps(&xre[i+48]); // a
                             zmm21 = _mm512_loadu_ps(&xim[i+48]); // b
                             const __mmask16 m3 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm18),
                                                              _mm512_abs_ps(zmm19),
                                                              _CMP_GE_OQ);
                             zmm22 = _mm512_mask_blend_ps(m3,_mm512_div_ps(zmm18,zmm19),
                                                             _mm512_div_ps(zmm19,zmm18)); // r = zmm22
                             zmm23 = _mm512_mask_blend_ps(m3,_mm512_fmadd_ps(zmm22,zmm18,zmm19),
                                                             _mm512_fmadd_ps(zmm22,zmm19,zmm18)); // den = zmm23
                             _mm512_storeu_ps(&zre[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm20,zmm22,zmm21),zmm23),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm21,zmm22,zmm20),zmm23)));
                             _mm512_storeu_ps(&zim[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm21,zmm22,zmm20),zmm23),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm21,_mm512_mul_ps(zmm20,zmm22)),zmm23)));
                             zmm24 = _mm512_loadu_ps(&yre[i+64]); // c
                             zmm25 = _mm512_loadu_ps(&yim[i+64]); // d
                             zmm26 = _mm512_loadu_ps(&xre[i+64]); // a
                             zmm27 = _mm512_loadu_ps(&xim[i+64]); // b
                             const __mmask16 m4 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm24),
                                                              _mm512_abs_ps(zmm25),
                                                              _CMP_GE_OQ);
                             zmm28 = _mm512_mask_blend_ps(m4,_mm512_div_ps(zmm24,zmm25),
                                                             _mm512_div_ps(zmm25,zmm24)); // r = zmm28
                             zmm29 = _mm512_mask_blend_ps(m4,_mm512_fmadd_ps(zmm28,zmm24,zmm25),
                                                             _mm512_fmadd_ps(zmm28,zmm25,zmm24)); // den = zmm29
                             _mm512_storeu_ps(&zre[i+64], _mm512_mask_blend_ps(m4,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm26,zmm28,zmm27),zmm29),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm27,zmm28,zmm26),zmm29)));
                             _mm512_storeu_ps(&zim[i+64], _mm512_mask_blend_ps(m4,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm27,zmm28,zmm26),zmm29),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm27,_mm512_mul_ps(zmm26,zmm28)),zmm29)));
                             zmm30 = _mm512_loadu_ps(&yre[i+80]); // c
                             zmm31 = _mm512_loadu_ps(&yim[i+80]); // d
                             zmm0  = _mm512_loadu_ps(&xre[i+80]); // a
                             zmm1  = _mm512_loadu_ps(&xim[i+80]); // b
                             const __mmask16 m5 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm30),
                                                              _mm512_abs_ps(zmm31),
                                                              _CMP_GE_OQ);
                             zmm2  = _mm512_mask_blend_ps(m5,_mm512_div_ps(zmm30,zmm31),
                                                             _mm512_div_ps(zmm31,zmm30)); // r = zmm2
                             zmm3  = _mm512_mask_blend_ps(m5,_mm512_fmadd_ps(zmm2,zmm30,zmm31),
                                                             _mm512_fmadd_ps(zmm2,zmm31,zmm30)); // den = zmm3
                             _mm512_storeu_ps(&zre[i+80], _mm512_mask_blend_ps(m5,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm0,zmm2,zmm1),zmm3),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm1,zmm2,zmm0),zmm3)));
                             _mm512_storeu_ps(&zim[i+80], _mm512_mask_blend_ps(m5,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm1,zmm2,zmm0),zmm3),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm1,_mm512_mul_ps(zmm0,zmm2)),zmm3)));
                        }

                          for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&yre[i+0]); // c
                             zmm1 = _mm512_loadu_ps(&yim[i+0]); // d
                             zmm2 = _mm512_loadu_ps(&xre[i+0]); // a
                             zmm3 = _mm512_loadu_ps(&xim[i+0]); // b
                             const __mmask16 m0 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                            _mm512_abs_ps(zmm1),
                                                            _CMP_GE_OQ);
                             zmm4 = _mm512_mask_blend_ps(m0,_mm512_div_ps(zmm0,zmm1),
                                                            _mm512_div_ps(zmm1,zmm0)); // r = zmm4
                             zmm5 = _mm512_mask_blend_ps(m0,_mm512_fmadd_ps(zmm4,zmm0,zmm1),
                                                            _mm512_fmadd_ps(zmm4,zmm1,zmm0)); // den = zmm5
                             _mm512_storeu_ps(&zre[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm2,zmm4,zmm3),zmm5),
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm3,zmm4,zmm2),zmm5)));
                             _mm512_storeu_ps(&zim[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmsub_ps(zmm3,zmm4,zmm2),zmm5),
                                                            _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,zmm4)),zmm5)));
                             zmm6 = _mm512_loadu_ps(&yre[i+16]); // c
                             zmm7 = _mm512_loadu_ps(&yim[i+16]); // d
                             zmm8 = _mm512_loadu_ps(&xre[i+16]); // a
                             zmm9 = _mm512_loadu_ps(&xim[i+16]); // b
                             const __mmask16 m1 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm6),
                                                              _mm512_abs_ps(zmm7),
                                                              _CMP_GE_OQ);
                             zmm10 = _mm512_mask_blend_ps(m1,_mm512_div_ps(zmm6,zmm7),
                                                             _mm512_div_ps(zmm7,zmm6)); // r = zmm10
                             zmm11 = _mm512_mask_blend_ps(m1,_mm512_fmadd_ps(zmm10,zmm6,zmm7),
                                                             _mm512_fmadd_ps(zmm10,zmm7,zmm6)); // den = zmm11
                             _mm512_storeu_ps(&zre[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm8,zmm10,zmm9),zmm11),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm9,zmm10,zmm8),zmm11)));
                             _mm512_storeu_ps(&zim[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm9,zmm10,zmm8),zmm11),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm9,_mm512_mul_ps(zmm8,zmm10)),zmm11)));
                             zmm12 = _mm512_loadu_ps(&yre[i+32]); // c
                             zmm13 = _mm512_loadu_ps(&yim[i+32]); // d
                             zmm14 = _mm512_loadu_ps(&xre[i+32]); // a
                             zmm15 = _mm512_loadu_ps(&xim[i+32]); // b
                             const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm12),
                                                              _mm512_abs_ps(zmm13),
                                                              _CMP_GE_OQ);
                             zmm16 = _mm512_mask_blend_ps(m2,_mm512_div_ps(zmm12,zmm13),
                                                             _mm512_div_ps(zmm13,zmm12)); // r = zmm16
                             zmm17 = _mm512_mask_blend_ps(m2,_mm512_fmadd_ps(zmm16,zmm12,zmm13),
                                                             _mm512_fmadd_ps(zmm16,zmm13,zmm12)); // den = zmm17
                             _mm512_storeu_ps(&zre[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm14,zmm16,zmm15),zmm17),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm15,zmm16,zmm14),zmm17)));
                             _mm512_storeu_ps(&zim[i+32], _mm512_mask_blend_ps(m2,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm15,zmm16,zmm14),zmm17),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm15,_mm512_mul_ps(zmm14,zmm16)),zmm17)));
                             zmm18 = _mm512_loadu_ps(&yre[i+48]); // c
                             zmm19 = _mm512_loadu_ps(&yim[i+48]); // d
                             zmm20 = _mm512_loadu_ps(&xre[i+48]); // a
                             zmm21 = _mm512_loadu_ps(&xim[i+48]); // b
                             const __mmask16 m3 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm18),
                                                              _mm512_abs_ps(zmm19),
                                                              _CMP_GE_OQ);
                             zmm22 = _mm512_mask_blend_ps(m3,_mm512_div_ps(zmm18,zmm19),
                                                             _mm512_div_ps(zmm19,zmm18)); // r = zmm22
                             zmm23 = _mm512_mask_blend_ps(m3,_mm512_fmadd_ps(zmm22,zmm18,zmm19),
                                                             _mm512_fmadd_ps(zmm22,zmm19,zmm18)); // den = zmm23
                             _mm512_storeu_ps(&zre[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm20,zmm22,zmm21),zmm23),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm21,zmm22,zmm20),zmm23)));
                             _mm512_storeu_ps(&zim[i+48], _mm512_mask_blend_ps(m3,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm21,zmm22,zmm20),zmm23),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm21,_mm512_mul_ps(zmm20,zmm22)),zmm23))); 
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&yre[i+0]); // c
                             zmm1 = _mm512_loadu_ps(&yim[i+0]); // d
                             zmm2 = _mm512_loadu_ps(&xre[i+0]); // a
                             zmm3 = _mm512_loadu_ps(&xim[i+0]); // b
                             const __mmask16 m0 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                            _mm512_abs_ps(zmm1),
                                                            _CMP_GE_OQ);
                             zmm4 = _mm512_mask_blend_ps(m0,_mm512_div_ps(zmm0,zmm1),
                                                            _mm512_div_ps(zmm1,zmm0)); // r = zmm4
                             zmm5 = _mm512_mask_blend_ps(m0,_mm512_fmadd_ps(zmm4,zmm0,zmm1),
                                                            _mm512_fmadd_ps(zmm4,zmm1,zmm0)); // den = zmm5
                             _mm512_storeu_ps(&zre[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm2,zmm4,zmm3),zmm5),
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm3,zmm4,zmm2),zmm5)));
                             _mm512_storeu_ps(&zim[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmsub_ps(zmm3,zmm4,zmm2),zmm5),
                                                            _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,zmm4)),zmm5)));
                             zmm6 = _mm512_loadu_ps(&yre[i+16]); // c
                             zmm7 = _mm512_loadu_ps(&yim[i+16]); // d
                             zmm8 = _mm512_loadu_ps(&xre[i+16]); // a
                             zmm9 = _mm512_loadu_ps(&xim[i+16]); // b
                             const __mmask16 m1 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm6),
                                                              _mm512_abs_ps(zmm7),
                                                              _CMP_GE_OQ);
                             zmm10 = _mm512_mask_blend_ps(m1,_mm512_div_ps(zmm6,zmm7),
                                                             _mm512_div_ps(zmm7,zmm6)); // r = zmm10
                             zmm11 = _mm512_mask_blend_ps(m1,_mm512_fmadd_ps(zmm10,zmm6,zmm7),
                                                             _mm512_fmadd_ps(zmm10,zmm7,zmm6)); // den = zmm11
                             _mm512_storeu_ps(&zre[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm8,zmm10,zmm9),zmm11),
                                                             _mm512_div_ps(_mm512_fmadd_ps(zmm9,zmm10,zmm8),zmm11)));
                             _mm512_storeu_ps(&zim[i+16], _mm512_mask_blend_ps(m1,
                                                             _mm512_div_ps(_mm512_fmsub_ps(zmm9,zmm10,zmm8),zmm11),
                                                             _mm512_div_ps(_mm512_sub_ps(zmm9,_mm512_mul_ps(zmm8,zmm10)),zmm11))); 
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&yre[i+0]); // c
                             zmm1 = _mm512_loadu_ps(&yim[i+0]); // d
                             zmm2 = _mm512_loadu_ps(&xre[i+0]); // a
                             zmm3 = _mm512_loadu_ps(&xim[i+0]); // b
                             const __mmask16 m0 = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                            _mm512_abs_ps(zmm1),
                                                            _CMP_GE_OQ);
                             zmm4 = _mm512_mask_blend_ps(m0,_mm512_div_ps(zmm0,zmm1),
                                                            _mm512_div_ps(zmm1,zmm0)); // r = zmm4
                             zmm5 = _mm512_mask_blend_ps(m0,_mm512_fmadd_ps(zmm4,zmm0,zmm1),
                                                            _mm512_fmadd_ps(zmm4,zmm1,zmm0)); // den = zmm5
                             _mm512_storeu_ps(&zre[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm2,zmm4,zmm3),zmm5),
                                                            _mm512_div_ps(_mm512_fmadd_ps(zmm3,zmm4,zmm2),zmm5)));
                             _mm512_storeu_ps(&zim[i+0], _mm512_mask_blend_ps(m0,
                                                            _mm512_div_ps(_mm512_fmsub_ps(zmm3,zmm4,zmm2),zmm5),
                                                            _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,zmm4)),zmm5)));
                        }

                          for(; (i+0) < n ; i += 1) {
                              const float yr = yre[i];
                              const float yi = yim[i];
                              const float xr = xre[i];
                              const float xi = xim[i];
                              
                              if(cephes_fabsf(yi) >= 
                                 cephes_fabsf(yr)) {
                                 float r   = yi/yr;
                                 float den = yr+r*yi;
                                 zre[i]    = (xr+xi*r)/den;
                                 zim[i]    = (xi-xr*r)/den;
                              }
                               else {
                                 float r   = yr/yi;
                                 float den = yi+r*yr;
                                 zre[i]    = (xr*r+xi)/den;
                                 zim[i]    = (xi*r+xr)/den;
                              }
                        }
 
               }


               /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





      } // math








} // gms





































#endif /*__GMS_CDIV_SMITH_VEC_ZMM16R4_HPP__*/
