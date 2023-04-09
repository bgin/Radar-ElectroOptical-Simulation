

#ifndef __GMS_CNORM_PROD_ZMM16R4_HPP__
#define __GMS_CNORM_PROD_ZMM16R4_HPP__


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

    const unsigned int GMS_CNORM_PROD_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CNORM_PROD_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CNORM_PROD_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CNORM_PROD_ZMM16R4_FULLVER =
      1000U*GMS_CNORM_PROD_ZMM16R4_MAJOR+
      100U*GMS_CNORM_PROD_ZMM16R4_MINOR+
      10U*GMS_CNORM_PROD_ZMM16R4_MICRO;
    const char * const GMS_CNORM_PROD_ZMM16R4_CREATION_DATE = "09-04-2023 09:44 AM +00200 ( SUN 09 APR 2023 GMT+2)";
    const char * const GMS_CNORM_PROD_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CNORM_PROD_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CNORM_PROD_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex norm product operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_cephes.h"


namespace  gms  {


          namespace  math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void cnorm_prod_u10x_zmm16r4_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  float * __restrict zre,
                                                  float * __restrict zim,
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
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_storeu_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_storeu_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_loadu_ps(&xre[i+80]);
                             zmm28 = _mm512_loadu_ps(&yre[i+80]);
                             zmm29 = _mm512_loadu_ps(&xim[i+80]);
                             zmm30 = _mm512_loadu_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm27 = _mm512_loadu_ps(&xre[i+96]);
                             zmm28 = _mm512_loadu_ps(&yre[i+96]);
                             zmm29 = _mm512_loadu_ps(&xim[i+96]);
                             zmm30 = _mm512_loadu_ps(&yim[i+96]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+96], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+96], _mm512_div_ps(zmm0,zmm29));
                             zmm1 = _mm512_loadu_ps(&xre[i+112]);
                             zmm2 = _mm512_loadu_ps(&yre[i+112]);
                             zmm3 = _mm512_loadu_ps(&xim[i+112]);
                             zmm4 = _mm512_loadu_ps(&yim[i+112]);
                             zmm5 = _mm512_fmsub_ps(zmm1,zmm2,
                                               _mm512_mul_ps(zmm3,zmm4)); // rep
                             zmm6  = _mm512_fmadd_pd(zmm3,zmm2,
                                               _mm512_mul_ps(zmm1,zmm4)); // imp
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_mul_ps(zmm6,zmm6);
                             zmm3 = _mm512_sqrt_ps(_mm512_add_ps(zmm5,zmm6)); // mag
                             _mm512_storeu_ps(&zre[i+112], _mm512_div_ps(zmm5,zmm3));
                             _mm512_storeu_ps(&zim[i+112], _mm512_div_ps(zmm6,zmm3));
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm4 = _mm512_loadu_ps(&xre[i+128]);
                             zmm5 = _mm512_loadu_ps(&yre[i+128]);
                             zmm6 = _mm512_loadu_ps(&xim[i+128]);
                             zmm7 = _mm512_loadu_ps(&yim[i+128]);
                             zmm8 = _mm512_fmsub_ps(zmm4,zmm5,
                                               _mm512_mul_ps(zmm6,zmm7)); // rep
                             zmm9  = _mm512_fmadd_pd(zmm6,zmm5,
                                               _mm512_mul_ps(zmm4,zmm7)); // imp
                             zmm4 = _mm512_mul_ps(zmm8,zmm8);
                             zmm5 = _mm512_mul_ps(zmm9,zmm9);
                             zmm6 = _mm512_sqrt_ps(_mm512_add_ps(zmm8,zmm9)); // mag
                             _mm512_storeu_ps(&zre[i+128], _mm512_div_ps(zmm8,zmm6));
                             _mm512_storeu_ps(&zim[i+128], _mm512_div_ps(zmm9,zmm6));
                             zmm10 = _mm512_loadu_ps(&xre[i+144]);
                             zmm11 = _mm512_loadu_ps(&yre[i+144]);
                             zmm12 = _mm512_loadu_ps(&xim[i+144]);
                             zmm13 = _mm512_loadu_ps(&yim[i+144]);
                             zmm14 = _mm512_fmsub_ps(zmm10,zmm11,
                                               _mm512_mul_ps(zmm12,zmm13)); // rep
                             zmm15  = _mm512_fmadd_pd(zmm12,zmm11,
                                               _mm512_mul_ps(zmm10,zmm13)); // imp
                             zmm10 = _mm512_mul_ps(zmm14,zmm14);
                             zmm11 = _mm512_mul_ps(zmm15,zmm15);
                             zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm10,zmm11)); // mag
                             _mm512_storeu_ps(&zre[i+144], _mm512_div_ps(zmm14,zmm12));
                             _mm512_storeu_ps(&zim[i+144], _mm512_div_ps(zmm15,zmm12));
                        }

                         for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_storeu_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_storeu_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_loadu_ps(&xre[i+80]);
                             zmm28 = _mm512_loadu_ps(&yre[i+80]);
                             zmm29 = _mm512_loadu_ps(&xim[i+80]);
                             zmm30 = _mm512_loadu_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                             zmm27 = _mm512_loadu_ps(&xre[i+96]);
                             zmm28 = _mm512_loadu_ps(&yre[i+96]);
                             zmm29 = _mm512_loadu_ps(&xim[i+96]);
                             zmm30 = _mm512_loadu_ps(&yim[i+96]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+96], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+96], _mm512_div_ps(zmm0,zmm29));
                             zmm1 = _mm512_loadu_ps(&xre[i+112]);
                             zmm2 = _mm512_loadu_ps(&yre[i+112]);
                             zmm3 = _mm512_loadu_ps(&xim[i+112]);
                             zmm4 = _mm512_loadu_ps(&yim[i+112]);
                             zmm5 = _mm512_fmsub_ps(zmm1,zmm2,
                                               _mm512_mul_ps(zmm3,zmm4)); // rep
                             zmm6  = _mm512_fmadd_pd(zmm3,zmm2,
                                               _mm512_mul_ps(zmm1,zmm4)); // imp
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_mul_ps(zmm6,zmm6);
                             zmm3 = _mm512_sqrt_ps(_mm512_add_ps(zmm5,zmm6)); // mag
                             _mm512_storeu_ps(&zre[i+112], _mm512_div_ps(zmm5,zmm3));
                             _mm512_storeu_ps(&zim[i+112], _mm512_div_ps(zmm6,zmm3));
                         
                        }

                          for(; (i+95) < n; i += 96) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_storeu_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_storeu_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_loadu_ps(&xre[i+80]);
                             zmm28 = _mm512_loadu_ps(&yre[i+80]);
                             zmm29 = _mm512_loadu_ps(&xim[i+80]);
                             zmm30 = _mm512_loadu_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                        }

                          for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));   
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                        }

                          for(; (i+1) < n; i += 1) {
                              const float xr  = xre[i];
                              const float yr  = yre[i];
                              const float xi  = xim[i];
                              const float yi  = yim[i];
                              const float re  = (xr*yr)-(xi*yi);
                              const float im  = (xi*yr)-(xr*yi);
                              const float mag = cephes_sqrtf(re*re+im*im);
                              zre[i]          = re/mag;
                              zim[i]          = im/mag;
                        }
                }


                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void cnorm_prod_u10x_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict __ATTR_ALIGN__(64) zre,
                                                  float * __restrict __ATTR_ALIGN__(64) zim,
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
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_store_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_store_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_store_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_store_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_load_ps(&xre[i+80]);
                             zmm28 = _mm512_load_ps(&yre[i+80]);
                             zmm29 = _mm512_load_ps(&xim[i+80]);
                             zmm30 = _mm512_load_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm27 = _mm512_load_ps(&xre[i+96]);
                             zmm28 = _mm512_load_ps(&yre[i+96]);
                             zmm29 = _mm512_load_ps(&xim[i+96]);
                             zmm30 = _mm512_load_ps(&yim[i+96]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+96], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+96], _mm512_div_ps(zmm0,zmm29));
                             zmm1 = _mm512_load_ps(&xre[i+112]);
                             zmm2 = _mm512_load_ps(&yre[i+112]);
                             zmm3 = _mm512_load_ps(&xim[i+112]);
                             zmm4 = _mm512_load_ps(&yim[i+112]);
                             zmm5 = _mm512_fmsub_ps(zmm1,zmm2,
                                               _mm512_mul_ps(zmm3,zmm4)); // rep
                             zmm6  = _mm512_fmadd_pd(zmm3,zmm2,
                                               _mm512_mul_ps(zmm1,zmm4)); // imp
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_mul_ps(zmm6,zmm6);
                             zmm3 = _mm512_sqrt_ps(_mm512_add_ps(zmm5,zmm6)); // mag
                             _mm512_store_ps(&zre[i+112], _mm512_div_ps(zmm5,zmm3));
                             _mm512_store_ps(&zim[i+112], _mm512_div_ps(zmm6,zmm3));
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm4 = _mm512_load_ps(&xre[i+128]);
                             zmm5 = _mm512_load_ps(&yre[i+128]);
                             zmm6 = _mm512_load_ps(&xim[i+128]);
                             zmm7 = _mm512_load_ps(&yim[i+128]);
                             zmm8 = _mm512_fmsub_ps(zmm4,zmm5,
                                               _mm512_mul_ps(zmm6,zmm7)); // rep
                             zmm9  = _mm512_fmadd_pd(zmm6,zmm5,
                                               _mm512_mul_ps(zmm4,zmm7)); // imp
                             zmm4 = _mm512_mul_ps(zmm8,zmm8);
                             zmm5 = _mm512_mul_ps(zmm9,zmm9);
                             zmm6 = _mm512_sqrt_ps(_mm512_add_ps(zmm8,zmm9)); // mag
                             _mm512_store_ps(&zre[i+128], _mm512_div_ps(zmm8,zmm6));
                             _mm512_store_ps(&zim[i+128], _mm512_div_ps(zmm9,zmm6));
                             zmm10 = _mm512_load_ps(&xre[i+144]);
                             zmm11 = _mm512_load_ps(&yre[i+144]);
                             zmm12 = _mm512_load_ps(&xim[i+144]);
                             zmm13 = _mm512_load_ps(&yim[i+144]);
                             zmm14 = _mm512_fmsub_ps(zmm10,zmm11,
                                               _mm512_mul_ps(zmm12,zmm13)); // rep
                             zmm15  = _mm512_fmadd_pd(zmm12,zmm11,
                                               _mm512_mul_ps(zmm10,zmm13)); // imp
                             zmm10 = _mm512_mul_ps(zmm14,zmm14);
                             zmm11 = _mm512_mul_ps(zmm15,zmm15);
                             zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm10,zmm11)); // mag
                             _mm512_store_ps(&zre[i+144], _mm512_div_ps(zmm14,zmm12));
                             _mm512_store_ps(&zim[i+144], _mm512_div_ps(zmm15,zmm12));
                        }

                         for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_store_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_store_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_load_ps(&xre[i+80]);
                             zmm28 = _mm512_load_ps(&yre[i+80]);
                             zmm29 = _mm512_load_ps(&xim[i+80]);
                             zmm30 = _mm512_load_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                             zmm27 = _mm512_load_ps(&xre[i+96]);
                             zmm28 = _mm512_load_ps(&yre[i+96]);
                             zmm29 = _mm512_load_ps(&xim[i+96]);
                             zmm30 = _mm512_load_ps(&yim[i+96]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+96], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+96], _mm512_div_ps(zmm0,zmm29));
                             zmm1 = _mm512_load_ps(&xre[i+112]);
                             zmm2 = _mm512_load_ps(&yre[i+112]);
                             zmm3 = _mm512_load_ps(&xim[i+112]);
                             zmm4 = _mm512_load_ps(&yim[i+112]);
                             zmm5 = _mm512_fmsub_ps(zmm1,zmm2,
                                               _mm512_mul_ps(zmm3,zmm4)); // rep
                             zmm6  = _mm512_fmadd_pd(zmm3,zmm2,
                                               _mm512_mul_ps(zmm1,zmm4)); // imp
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_mul_ps(zmm6,zmm6);
                             zmm3 = _mm512_sqrt_ps(_mm512_add_ps(zmm5,zmm6)); // mag
                             _mm512_store_ps(&zre[i+112], _mm512_div_ps(zmm5,zmm3));
                             _mm512_store_ps(&zim[i+112], _mm512_div_ps(zmm6,zmm3));
                         
                        }

                          for(; (i+95) < n; i += 96) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_store_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_store_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_store_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_store_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_load_ps(&xre[i+80]);
                             zmm28 = _mm512_load_ps(&yre[i+80]);
                             zmm29 = _mm512_load_ps(&xim[i+80]);
                             zmm30 = _mm512_load_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                        }

                          for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_store_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_store_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));   
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                        }

                          for(; (i+1) < n; i += 1) {
                              const float xr  = xre[i];
                              const float yr  = yre[i];
                              const float xi  = xim[i];
                              const float yi  = yim[i];
                              const float re  = (xr*yr)-(xi*yi);
                              const float im  = (xi*yr)-(xr*yi);
                              const float mag = cephes_sqrtf(re*re+im*im);
                              zre[i]          = re/mag;
                              zim[i]          = im/mag;
                        }
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void cnorm_prod_u8x_zmm16r4_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  float * __restrict zre,
                                                  float * __restrict zim,
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
                         
                         for(i = 0; (i+127) < n; i += 128) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_storeu_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_storeu_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_loadu_ps(&xre[i+80]);
                             zmm28 = _mm512_loadu_ps(&yre[i+80]);
                             zmm29 = _mm512_loadu_ps(&xim[i+80]);
                             zmm30 = _mm512_loadu_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm27 = _mm512_loadu_ps(&xre[i+96]);
                             zmm28 = _mm512_loadu_ps(&yre[i+96]);
                             zmm29 = _mm512_loadu_ps(&xim[i+96]);
                             zmm30 = _mm512_loadu_ps(&yim[i+96]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+96], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+96], _mm512_div_ps(zmm0,zmm29));
                             zmm1 = _mm512_loadu_ps(&xre[i+112]);
                             zmm2 = _mm512_loadu_ps(&yre[i+112]);
                             zmm3 = _mm512_loadu_ps(&xim[i+112]);
                             zmm4 = _mm512_loadu_ps(&yim[i+112]);
                             zmm5 = _mm512_fmsub_ps(zmm1,zmm2,
                                               _mm512_mul_ps(zmm3,zmm4)); // rep
                             zmm6  = _mm512_fmadd_pd(zmm3,zmm2,
                                               _mm512_mul_ps(zmm1,zmm4)); // imp
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_mul_ps(zmm6,zmm6);
                             zmm3 = _mm512_sqrt_ps(_mm512_add_ps(zmm5,zmm6)); // mag
                             _mm512_storeu_ps(&zre[i+112], _mm512_div_ps(zmm5,zmm3));
                             _mm512_storeu_ps(&zim[i+112], _mm512_div_ps(zmm6,zmm3));
                          
                        }

                          for(; (i+95) < n; i += 96) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_storeu_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_storeu_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_loadu_ps(&xre[i+80]);
                             zmm28 = _mm512_loadu_ps(&yre[i+80]);
                             zmm29 = _mm512_loadu_ps(&xim[i+80]);
                             zmm30 = _mm512_loadu_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                        }

                          for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));   
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                        }

                          for(; (i+1) < n; i += 1) {
                              const float xr  = xre[i];
                              const float yr  = yre[i];
                              const float xi  = xim[i];
                              const float yi  = yim[i];
                              const float re  = (xr*yr)-(xi*yi);
                              const float im  = (xi*yr)-(xr*yi);
                              const float mag = cephes_sqrtf(re*re+im*im);
                              zre[i]          = re/mag;
                              zim[i]          = im/mag;
                        }
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void cnorm_prod_u8x_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict __ATTR_ALIGN__(64) zre,
                                                  float * __restrict __ATTR_ALIGN__(64) zim,
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
                         
                         for(i = 0; (i+127) < n; i += 128) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_store_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_store_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_store_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_store_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_load_ps(&xre[i+80]);
                             zmm28 = _mm512_load_ps(&yre[i+80]);
                             zmm29 = _mm512_load_ps(&xim[i+80]);
                             zmm30 = _mm512_load_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm27 = _mm512_load_ps(&xre[i+96]);
                             zmm28 = _mm512_load_ps(&yre[i+96]);
                             zmm29 = _mm512_load_ps(&xim[i+96]);
                             zmm30 = _mm512_load_ps(&yim[i+96]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+96], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+96], _mm512_div_ps(zmm0,zmm29));
                             zmm1 = _mm512_load_ps(&xre[i+112]);
                             zmm2 = _mm512_load_ps(&yre[i+112]);
                             zmm3 = _mm512_load_ps(&xim[i+112]);
                             zmm4 = _mm512_load_ps(&yim[i+112]);
                             zmm5 = _mm512_fmsub_ps(zmm1,zmm2,
                                               _mm512_mul_ps(zmm3,zmm4)); // rep
                             zmm6  = _mm512_fmadd_pd(zmm3,zmm2,
                                               _mm512_mul_ps(zmm1,zmm4)); // imp
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_mul_ps(zmm6,zmm6);
                             zmm3 = _mm512_sqrt_ps(_mm512_add_ps(zmm5,zmm6)); // mag
                             _mm512_store_ps(&zre[i+112], _mm512_div_ps(zmm5,zmm3));
                             _mm512_store_ps(&zim[i+112], _mm512_div_ps(zmm6,zmm3));
                        }

                         for(; (i+95) < n; i += 96) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_store_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_store_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_store_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_store_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_load_ps(&xre[i+80]);
                             zmm28 = _mm512_load_ps(&yre[i+80]);
                             zmm29 = _mm512_load_ps(&xim[i+80]);
                             zmm30 = _mm512_load_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_store_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_store_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                        }

                          for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_load_ps(&xre[i+32]);
                             zmm13 = _mm512_load_ps(&yre[i+32]);
                             zmm14 = _mm512_load_ps(&xim[i+32]);
                             zmm15 = _mm512_load_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_store_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_store_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_store_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_store_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));   
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_load_ps(&xre[i+16]);
                             zmm7 = _mm512_load_ps(&yre[i+16]);
                             zmm8 = _mm512_load_ps(&xim[i+16]);
                             zmm9 = _mm512_load_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_store_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_store_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_load_ps(&xre[i+0]);
                             zmm1 = _mm512_load_ps(&yre[i+0]);
                             zmm2 = _mm512_load_ps(&xim[i+0]);
                             zmm3 = _mm512_load_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_store_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_store_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                        }

                          for(; (i+1) < n; i += 1) {
                              const float xr  = xre[i];
                              const float yr  = yre[i];
                              const float xi  = xim[i];
                              const float yi  = yim[i];
                              const float re  = (xr*yr)-(xi*yi);
                              const float im  = (xi*yr)-(xr*yi);
                              const float mag = cephes_sqrtf(re*re+im*im);
                              zre[i]          = re/mag;
                              zim[i]          = im/mag;
                        }
                }


                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void cnorm_prod_u6x_zmm16r4_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  float * __restrict zre,
                                                  float * __restrict zim,
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
                         
                         for(i = 0; (i+95) < n; i += 96) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             zmm29 = _mm512_fmadd_pd(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             zmm24 = _mm512_mul_ps(zmm28,zmm28);
                             zmm25 = _mm512_mul_ps(zmm29,zmm29);
                             zmm26 = _mm512_sqrt_ps(_mm512_add_ps(zmm24,zmm25)); // mag
                             _mm512_storeu_ps(&zre[i+64], _mm512_div_ps(zmm28,zmm26));
                             _mm512_storeu_ps(&zim[i+64], _mm512_div_ps(zmm29,zmm26));
                             zmm27 = _mm512_loadu_ps(&xre[i+80]);
                             zmm28 = _mm512_loadu_ps(&yre[i+80]);
                             zmm29 = _mm512_loadu_ps(&xim[i+80]);
                             zmm30 = _mm512_loadu_ps(&yim[i+80]);
                             zmm31 = _mm512_fmsub_ps(zmm27,zmm28,
                                               _mm512_mul_ps(zmm29,zmm30)); // rep
                             zmm0  = _mm512_fmadd_pd(zmm29,zmm28,
                                               _mm512_mul_ps(zmm27,zmm30)); // imp
                             zmm27 = _mm512_mul_ps(zmm31,zmm31);
                             zmm28 = _mm512_mul_ps(zmm0,zmm0);
                             zmm29 = _mm512_sqrt_ps(_mm512_add_ps(zmm27,zmm28)); // mag
                             _mm512_storeu_ps(&zre[i+80], _mm512_div_ps(zmm31,zmm29));
                             _mm512_storeu_ps(&zim[i+80], _mm512_div_ps(zmm0,zmm29));
                                                      
                        }

                        for(; (i+63) < n; i += 64) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                             zmm12 = _mm512_loadu_ps(&xre[i+32]);
                             zmm13 = _mm512_loadu_ps(&yre[i+32]);
                             zmm14 = _mm512_loadu_ps(&xim[i+32]);
                             zmm15 = _mm512_loadu_ps(&yim[i+32]);
                             zmm16= _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             zmm17= _mm512_fmadd_pd(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             zmm12 = _mm512_mul_ps(zmm16,zmm16);
                             zmm13 = _mm512_mul_ps(zmm17,zmm17);
                             zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm12,zmm13));
                             _mm512_storeu_ps(&zre[i+32], _mm512_div_ps(zmm16,zmm14));
                             _mm512_storeu_ps(&zim[i+32], _mm512_div_ps(zmm17,zmm14));
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             zmm23 = _mm512_fmadd_pd(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             zmm18 = _mm512_mul_ps(zmm22,zmm22);
                             zmm19 = _mm512_mul_ps(zmm23,zmm23);
                             zmm20 = _mm512_sqrt_ps(_mm512_add_ps(zmm18,zmm19)); // mag
                             _mm512_storeu_ps(&zre[i+48], _mm512_div_ps(zmm22,zmm20));
                             _mm512_storeu_ps(&zim[i+48], _mm512_div_ps(zmm23,zmm20));   
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                             zmm6 = _mm512_loadu_ps(&xre[i+16]);
                             zmm7 = _mm512_loadu_ps(&yre[i+16]);
                             zmm8 = _mm512_loadu_ps(&xim[i+16]);
                             zmm9 = _mm512_loadu_ps(&yim[i+16]);
                             zmm10= _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             zmm11= _mm512_fmadd_pd(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             zmm6 = _mm512_mul_ps(zmm10,zmm10);
                             zmm7 = _mm512_mul_ps(zmm11,zmm11);
                             zmm8 = _mm512_sqrt_ps(_mm512_add_ps(zmm6,zmm7));
                             _mm512_storeu_ps(&zre[i+16], _mm512_div_ps(zmm10,zmm8));
                             _mm512_storeu_ps(&zim[i+16], _mm512_div_ps(zmm11,zmm8));
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0 = _mm512_loadu_ps(&xre[i+0]);
                             zmm1 = _mm512_loadu_ps(&yre[i+0]);
                             zmm2 = _mm512_loadu_ps(&xim[i+0]);
                             zmm3 = _mm512_loadu_ps(&yim[i+0]);
                             zmm4 = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             zmm5 = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             zmm0 = _mm512_mul_ps(zmm4,zmm4);
                             zmm1 = _mm512_mul_ps(zmm5,zmm5);
                             zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                             _mm512_storeu_ps(&zre[i+0], _mm512_div_ps(zmm4,zmm2));
                             _mm512_storeu_ps(&zim[i+0], _mm512_div_ps(zmm5,zmm2));
                        }

                          for(; (i+1) < n; i += 1) {
                              const float xr  = xre[i];
                              const float yr  = yre[i];
                              const float xi  = xim[i];
                              const float yi  = yim[i];
                              const float re  = (xr*yr)-(xi*yi);
                              const float im  = (xi*yr)-(xr*yi);
                              const float mag = cephes_sqrtf(re*re+im*im);
                              zre[i]          = re/mag;
                              zim[i]          = im/mag;
                        }
                }



       } // math


} // gms















#endif /*__GMS_CNORM_PROD_ZMM16R4_HPP__*/
