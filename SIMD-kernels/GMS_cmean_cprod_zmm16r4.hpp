

#ifndef __GMS_CMEAN_CPROD_ZMM16R4_HPP__
#define __GMS_CMEAN_CPROD_ZMM16R4_HPP__ 100420230944

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

    const unsigned int GMS_CMEAN_CPROD_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CMEAN_CPROD_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CMEAN_CPROD_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CMEAN_CPROD_ZMM16R4_FULLVER =
      1000U*GMS_CMEAN_CPROD_ZMM16R4_MAJOR+
      100U*GMS_CMEAN_CPROD_ZMM16R4_MINOR+
      10U*GMS_CMEAN_CPROD_ZMM16R4_MICRO;
    const char * const GMS_CMEAN_CPROD_ZMM16R4_CREATION_DATE = "09-04-2023 09:44 AM +00200 ( SUN 09 APR 2023 GMT+2)";
    const char * const GMS_CMEAN_CPROD_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CMEAN_CPROD_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CMEAN_CPROD_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex mean product operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_zmm16r4.hpp"


namespace gms {


         namespace math {


                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline     
                   void cmean_prod_u10x_zmm16r4_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  float * __restrict mre,
                                                  float * __restrict mim,
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
                         register float sre,sim,accr,acci;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         sre = 0.0f;
                         sim = 0.0f;
                         accr= 0.0f;
                         acci= 0.0f;
                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             sre   = _mm512_reduce_add_ps(zmm10);
                             accr += sre;
                             zmm11 = _mm512_fmadd_ps(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             sim   = _mm512_reduce_add_ps(zmm11);
                             acci  += sim;
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             sre   = _mm512_reduce_add_ps(zmm16);
                             accr += sre;
                             zmm17 = _mm512_fmadd_ps(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             sim   = _mm512_reduce_add_ps(zmm17);
                             acci  += sim;
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             sre   = _mm512_reduce_add_ps(zmm22);
                             accr += sre;
                             zmm23 = _mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             sim   = _mm512_reduce_add_ps(zmm23);
                             acci  += sim;
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
                             sre   = _mm512_reduce_add_ps(zmm28);
                             accr += sre;
                             zmm29 = _mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             sim   = _mm512_reduce_add_ps(zmm29);
                             acci  += sim;
                             zmm30 = _mm512_loadu_ps(&xre[i+80]);
                             zmm31 = _mm512_loadu_ps(&yre[i+80]);
                             zmm0  = _mm512_loadu_ps(&xim[i+80]);
                             zmm1  = _mm512_loadu_ps(&yim[i+80]);
                             zmm2  = _mm512_fmsub_ps(zmm30,zmm31,
                                               _mm512_mul_ps(zmm0,zmm1)); // rep
                             sre   = _mm512_reduce_add_ps(zmm2);
                             accr  += sre;
                             zmm3  = _mm512_fmadd_ps(zmm0,zmm31,
                                               _mm512_mul_ps(zmm30,zmm1)); // imp
                             sim   = _mm512_reduce_add_ps(zmm3);
                             acci  += sim;
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm4  = _mm512_loadu_ps(&xre[i+96]);
                             zmm5  = _mm512_loadu_ps(&yre[i+96]);
                             zmm6  = _mm512_loadu_ps(&xim[i+96]);
                             zmm7  = _mm512_loadu_ps(&yim[i+96]);
                             zmm8  = _mm512_fmsub_ps(zmm4,zmm5,
                                               _mm512_mul_ps(zmm6,zmm7)); // rep
                             sre   = _mm512_reduce_add_ps(zmm8);
                             accr  += sre;
                             zmm9  = _mm512_fmadd_ps(zmm6,zmm5,
                                               _mm512_mul_ps(zmm4,zmm7)); // imp
                             sim   = _mm512_reduce_add_ps(zmm9);
                             acci  += sim;
                             zmm10 = _mm512_loadu_ps(&xre[i+112]);
                             zmm11 = _mm512_loadu_ps(&yre[i+112]);
                             zmm12 = _mm512_loadu_ps(&xim[i+112]);
                             zmm13 = _mm512_loadu_ps(&yim[i+112]);
                             zmm14 = _mm512_fmsub_ps(zmm10,zmm11,
                                               _mm512_mul_ps(zmm12,zmm13)); // rep
                             sre   = _mm512_reduce_add_ps(zmm14);
                             accr  += sre;
                             zmm15 = _mm512_fmadd_ps(zmm13,zmm12,
                                               _mm512_mul_ps(zmm10,zmm14)); // imp
                             sim   = _mm512_reduce_add_ps(zmm15);
                             acci  += sim;
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm16 = _mm512_loadu_ps(&xre[i+128]);
                             zmm17 = _mm512_loadu_ps(&yre[i+128]);
                             zmm18 = _mm512_loadu_ps(&xim[i+128]);
                             zmm19 = _mm512_loadu_ps(&yim[i+128]);
                             zmm20 = _mm512_fmsub_ps(zmm16,zmm17,
                                               _mm512_mul_ps(zmm18,zmm19)); // rep
                             sre   = _mm512_reduce_add_ps(zmm20);
                             accr  += sre;
                             zmm21 = _mm512_fmadd_ps(zmm18,zmm17,
                                               _mm512_mul_ps(zmm16,zmm19)); // imp
                             sim   = _mm512_reduce_add_ps(zmm21);
                             acci  += sim;
                             zmm22 = _mm512_loadu_ps(&xre[i+144]);
                             zmm23 = _mm512_loadu_ps(&yre[i+144]);
                             zmm24 = _mm512_loadu_ps(&xim[i+144]);
                             zmm25 = _mm512_loadu_ps(&yim[i+144]);
                             zmm26 = _mm512_fmsub_ps(zmm22,zmm23,
                                               _mm512_mul_ps(zmm24,zmm25)); // rep
                             sre   = _mm512_reduce_add_ps(zmm26);
                             accr  += sre;
                             zmm27 = _mm512_fmadd_ps(zmm24,zmm23,
                                               _mm512_mul_ps(zmm22,zmm25)); // imp
                             sim   = _mm512_reduce_add_ps(zmm27);
                             acci  += sim;
                        }

                         for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             sre   = _mm512_reduce_add_ps(zmm10);
                             accr += sre;
                             zmm11 = _mm512_fmadd_ps(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             sim   = _mm512_reduce_add_ps(zmm11);
                             acci  += sim;
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             sre   = _mm512_reduce_add_ps(zmm16);
                             accr += sre;
                             zmm17 = _mm512_fmadd_ps(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             sim   = _mm512_reduce_add_ps(zmm17);
                             acci  += sim;
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             sre   = _mm512_reduce_add_ps(zmm22);
                             accr += sre;
                             zmm23 = _mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             sim   = _mm512_reduce_add_ps(zmm23);
                             acci  += sim;
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             sre   = _mm512_reduce_add_ps(zmm28);
                             accr += sre;
                             zmm29 = _mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             sim   = _mm512_reduce_add_ps(zmm29);
                             acci  += sim;
                             zmm30 = _mm512_loadu_ps(&xre[i+80]);
                             zmm31 = _mm512_loadu_ps(&yre[i+80]);
                             zmm0  = _mm512_loadu_ps(&xim[i+80]);
                             zmm1  = _mm512_loadu_ps(&yim[i+80]);
                             zmm2  = _mm512_fmsub_ps(zmm30,zmm31,
                                               _mm512_mul_ps(zmm0,zmm1)); // rep
                             sre   = _mm512_reduce_add_ps(zmm2);
                             accr  += sre;
                             zmm3  = _mm512_fmadd_ps(zmm0,zmm31,
                                               _mm512_mul_ps(zmm30,zmm1)); // imp
                             sim   = _mm512_reduce_add_ps(zmm3);
                             acci  += sim;
                             zmm4  = _mm512_loadu_ps(&xre[i+96]);
                             zmm5  = _mm512_loadu_ps(&yre[i+96]);
                             zmm6  = _mm512_loadu_ps(&xim[i+96]);
                             zmm7  = _mm512_loadu_ps(&yim[i+96]);
                             zmm8  = _mm512_fmsub_ps(zmm4,zmm5,
                                               _mm512_mul_ps(zmm6,zmm7)); // rep
                             sre   = _mm512_reduce_add_ps(zmm8);
                             accr  += sre;
                             zmm9  = _mm512_fmadd_ps(zmm6,zmm5,
                                               _mm512_mul_ps(zmm4,zmm7)); // imp
                             sim   = _mm512_reduce_add_ps(zmm9);
                             acci  += sim;
                             zmm10 = _mm512_loadu_ps(&xre[i+112]);
                             zmm11 = _mm512_loadu_ps(&yre[i+112]);
                             zmm12 = _mm512_loadu_ps(&xim[i+112]);
                             zmm13 = _mm512_loadu_ps(&yim[i+112]);
                             zmm14 = _mm512_fmsub_ps(zmm10,zmm11,
                                               _mm512_mul_ps(zmm12,zmm13)); // rep
                             sre   = _mm512_reduce_add_ps(zmm14);
                             accr  += sre;
                             zmm15 = _mm512_fmadd_ps(zmm13,zmm12,
                                               _mm512_mul_ps(zmm10,zmm14)); // imp
                             sim   = _mm512_reduce_add_ps(zmm15);
                             acci  += sim;
                        }

                          for(; (i+95) < n; i += 96) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             sre   = _mm512_reduce_add_ps(zmm10);
                             accr += sre;
                             zmm11 = _mm512_fmadd_ps(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             sim   = _mm512_reduce_add_ps(zmm11);
                             acci  += sim;
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             sre   = _mm512_reduce_add_ps(zmm16);
                             accr += sre;
                             zmm17 = _mm512_fmadd_ps(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             sim   = _mm512_reduce_add_ps(zmm17);
                             acci  += sim;
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             sre   = _mm512_reduce_add_ps(zmm22);
                             accr += sre;
                             zmm23 = _mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             sim   = _mm512_reduce_add_ps(zmm23);
                             acci  += sim;
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             sre   = _mm512_reduce_add_ps(zmm28);
                             accr += sre;
                             zmm29 = _mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             sim   = _mm512_reduce_add_ps(zmm29);
                             acci  += sim;
                             zmm30 = _mm512_loadu_ps(&xre[i+80]);
                             zmm31 = _mm512_loadu_ps(&yre[i+80]);
                             zmm0  = _mm512_loadu_ps(&xim[i+80]);
                             zmm1  = _mm512_loadu_ps(&yim[i+80]);
                             zmm2  = _mm512_fmsub_ps(zmm30,zmm31,
                                               _mm512_mul_ps(zmm0,zmm1)); // rep
                             sre   = _mm512_reduce_add_ps(zmm2);
                             accr  += sre;
                             zmm3  = _mm512_fmadd_ps(zmm0,zmm31,
                                               _mm512_mul_ps(zmm30,zmm1)); // imp
                             sim   = _mm512_reduce_add_ps(zmm3);
                             acci  += sim;
                        }

                          for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             sre   = _mm512_reduce_add_ps(zmm10);
                             accr += sre;
                             zmm11 = _mm512_fmadd_ps(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             sim   = _mm512_reduce_add_ps(zmm11);
                             acci  += sim;
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             sre   = _mm512_reduce_add_ps(zmm16);
                             accr += sre;
                             zmm17 = _mm512_fmadd_ps(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             sim   = _mm512_reduce_add_ps(zmm17);
                             acci  += sim;
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             sre   = _mm512_reduce_add_ps(zmm22);
                             accr += sre;
                             zmm23 = _mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             sim   = _mm512_reduce_add_ps(zmm23);
                             acci  += sim;
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             zmm28 = _mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27)); // rep
                             sre   = _mm512_reduce_add_ps(zmm28);
                             accr += sre;
                             zmm29 = _mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27)); // imp
                             sim   = _mm512_reduce_add_ps(zmm29);
                             acci  += sim;
                        }

                          for(; (i+63) < n; i += 64) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             sre   = _mm512_reduce_add_ps(zmm10);
                             accr += sre;
                             zmm11 = _mm512_fmadd_ps(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             sim   = _mm512_reduce_add_ps(zmm11);
                             acci  += sim;
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             zmm16 = _mm512_fmsub_ps(zmm12,zmm13,
                                               _mm512_mul_ps(zmm14,zmm15)); // rep
                             sre   = _mm512_reduce_add_ps(zmm16);
                             accr += sre;
                             zmm17 = _mm512_fmadd_ps(zmm14,zmm13,
                                               _mm512_mul_ps(zmm12,zmm15)); // imp
                             sim   = _mm512_reduce_add_ps(zmm17);
                             acci  += sim;
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             zmm22 = _mm512_fmsub_ps(zmm18,zmm19,
                                               _mm512_mul_ps(zmm20,zmm21)); // rep
                             sre   = _mm512_reduce_add_ps(zmm22);
                             accr += sre;
                             zmm23 = _mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21)); // imp
                             sim   = _mm512_reduce_add_ps(zmm23);
                             acci  += sim;
                        }

                          for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             zmm10 = _mm512_fmsub_ps(zmm6,zmm7,
                                               _mm512_mul_ps(zmm8,zmm9)); // rep
                             sre   = _mm512_reduce_add_ps(zmm10);
                             accr += sre;
                             zmm11 = _mm512_fmadd_ps(zmm8,zmm7,
                                               _mm512_mul_ps(zmm6,zmm9)); // imp
                             sim   = _mm512_reduce_add_ps(zmm11);
                             acci  += sim;
                        }

                          for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             zmm4  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3)); // rep
                             sre   = _mm512_reduce_add_ps(zmm4);
                             accr  += sre;
                             zmm5  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3)); // imp
                             sim   = _mm512_reduce_add_ps(zmm5);
                             acci  += sim;
                        }

                          for(; (i+0) < n; i += 1) {
                              const float xr = xre[i];
                              const float yr = yre[i];
                              const float xi = xim[i];
                              const float yi = yim[i];
                              const float re = (xr*yr)-(xi*yi);
                              const float im = (xi*yr)+(xr*yi);
                              accr += re;
                              acci += im;
                        }
                          *mre = accr*invN;
                          *mim = acci*invN;
                } 
               
      } // math

} // gms













#endif /*__GMS_CMEAN_CPROD_ZMM16R4*/
