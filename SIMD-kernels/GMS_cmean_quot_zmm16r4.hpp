

#ifndef __GMS_CMEAN_QUOT_ZMM16R4_HPP__
#define __GMS_CMEAN_QUOT_ZMM16R4_HPP__ 120120230939


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

    const unsigned int GMS_CMEAN_QUOT_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CMEAN_QUOT_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CMEAN_QUOT_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CMEAN_QUOT_ZMM16R4_FULLVER =
      1000U*GMS_CMEAN_QUOT_ZMM16R4_MAJOR+
      100U*GMS_CMEAN_QUOT_ZMM16R4_MINOR+
      10U*GMS_CMEAN_QUOT_ZMM16R4_MICRO;
    const char * const GMS_CMEAN_QUOT_ZMM16R4_CREATION_DATE = "12-04-2023 09:39 AM +00200 ( WED 12 APR 2023 GMT+2)";
    const char * const GMS_CMEAN_QUOT_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CMEAN_QUOT_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CMEAN_QUOT_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex mean product operations."

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
                   void cmean_quot_u10x_zmm16r4_u(const float * __restrict xre,
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
                         __m512 redr[10] = {_mm512_setzero_ps()};
                         __m512 redi[10] = {_mm512_setzero_ps()};
                         register float re,im;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         re = 0.0f;
                         im = 0.0f;

                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0    = _mm512_loadu_ps(&xre[i+0]);
                             zmm1    = _mm512_loadu_ps(&yre[i+0]);
                             zmm2    = _mm512_loadu_ps(&xim[i+0]);
                             zmm3    = _mm512_loadu_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4));
                             zmm7    = _mm512_loadu_ps(&xre[i+16]);
                             zmm8    = _mm512_loadu_ps(&yre[i+16]);
                             zmm9    = _mm512_loadu_ps(&xim[i+16]);
                             zmm10   = _mm512_loadu_ps(&yim[i+16]);
                             zmm11   = _mm512_fmadd_ps(zmm8,zmm8,
                                                          _mm512_mul_ps(zmm10,zmm10)); // den
                             zmm12   = _mm512_fmsub_ps(zmm7,zmm8,
                                                          _mm512_mul_ps(zmm9,zmm10));// rep
                             redr[1] = _mm512_add_ps(redr[1],_mm512_div_ps(zmm12,zmm11)); 
                             zmm13   = _mm512_fmadd_ps(zmm9,zmm8,
                                                          _mm512_mul_ps(zmm7,zmm10)); // imp
                             redi[1] = _mm512_add_ps(redi[1],_mm512_div_ps(zmm13,zmm11));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm14   = _mm512_loadu_ps(&xre[i+32]);
                             zmm15   = _mm512_loadu_ps(&yre[i+32]);
                             zmm16   = _mm512_loadu_ps(&xim[i+32]);
                             zmm17   = _mm512_loadu_ps(&yim[i+32]);
                             zmm18   = _mm512_fmadd_ps(zmm15,zmm15,
                                                          _mm512_mul_ps(zmm17,zmm17)); // den
                             zmm19   = _mm512_fmsub_ps(zmm14,zmm15,
                                                          _mm512_mul_ps(zmm16,zmm17));// rep
                             redr[2] = _mm512_add_ps(redr[2],_mm512_div_ps(zmm19,zmm18)); 
                             zmm20   = _mm512_fmadd_ps(zmm16,zmm15,
                                                          _mm512_mul_ps(zmm14,zmm17)); // imp
                             redi[2] = _mm512_add_ps(redi[2],_mm512_div_ps(zmm20,zmm18));  
                             zmm21   = _mm512_loadu_ps(&xre[i+48]);
                             zmm22   = _mm512_loadu_ps(&yre[i+48]);
                             zmm23   = _mm512_loadu_ps(&xim[i+48]);
                             zmm24   = _mm512_loadu_ps(&yim[i+48]);
                             zmm25   = _mm512_fmadd_ps(zmm22,zmm22,
                                                          _mm512_mul_ps(zmm24,zmm24)); // den
                             zmm26   = _mm512_fmsub_ps(zmm21,zmm22,
                                                          _mm512_mul_ps(zmm23,zmm24));// rep
                             redr[3] = _mm512_add_ps(redr[3],_mm512_div_ps(zmm26,zmm25)); 
                             zmm27   = _mm512_fmadd_ps(zmm23,zmm22,
                                                          _mm512_mul_ps(zmm21,zmm24)); // imp
                             redi[3] = _mm512_add_ps(redi[3],_mm512_div_ps(zmm27,zmm25)); 
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0); 
                             zmm28   = _mm512_loadu_ps(&xre[i+64]);
                             zmm29   = _mm512_loadu_ps(&yre[i+64]);
                             zmm30   = _mm512_loadu_ps(&xim[i+64]);
                             zmm31   = _mm512_loadu_ps(&yim[i+64]);
                             zmm0    = _mm512_fmadd_ps(zmm29,zmm29,
                                                          _mm512_mul_ps(zmm31,zmm31)); // den
                             zmm1    = _mm512_fmsub_ps(zmm28,zmm29,
                                                          _mm512_mul_ps(zmm30,zmm31));// rep
                             redr[4] = _mm512_add_ps(redr[4],_mm512_div_ps(zmm1,zmm0)); 
                             zmm2    = _mm512_fmadd_ps(zmm30,zmm29,
                                                          _mm512_mul_ps(zmm28,zmm31)); // imp
                             redi[4] = _mm512_add_ps(redi[4],_mm512_div_ps(zmm2,zmm0));   
                             zmm3   = _mm512_loadu_ps(&xre[i+80]);
                             zmm4   = _mm512_loadu_ps(&yre[i+80]);
                             zmm5   = _mm512_loadu_ps(&xim[i+80]);
                             zmm6   = _mm512_loadu_ps(&yim[i+80]);
                             zmm7   = _mm512_fmadd_ps(zmm4,zmm4,
                                                          _mm512_mul_ps(zmm6,zmm6)); // den
                             zmm8   = _mm512_fmsub_ps(zmm3,zmm4,
                                                          _mm512_mul_ps(zmm5,zmm6));// rep
                             redr[5]= _mm512_add_ps(redr[5],_mm512_div_ps(zmm8,zmm7)); 
                             zmm9   = _mm512_fmadd_ps(zmm5,zmm4,
                                                          _mm512_mul_ps(zmm3,zmm6)); // imp
                             redi[5]= _mm512_add_ps(redi[5],_mm512_div_ps(zmm9,zmm7)); 
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm10  = _mm512_loadu_ps(&xre[i+96]);
                             zmm11  = _mm512_loadu_ps(&yre[i+96]);
                             zmm12  = _mm512_loadu_ps(&xim[i+96]);
                             zmm13  = _mm512_loadu_ps(&yim[i+96]);
                             zmm14  = _mm512_fmadd_ps(zmm11,zmm11,
                                                          _mm512_mul_ps(zmm13,zmm13)); // den
                             zmm15   = _mm512_fmsub_ps(zmm10,zmm11,
                                                          _mm512_mul_ps(zmm12,zmm13));// rep
                             redr[6] = _mm512_add_ps(redr[6],_mm512_div_ps(zmm15,zmm14)); 
                             zmm16    = _mm512_fmadd_ps(zmm12,zmm11,
                                                          _mm512_mul_ps(zmm10,zmm13)); // imp
                             redi[6] = _mm512_add_ps(redi[6],_mm512_div_ps(zmm16,zmm14)); 
                             zmm17  = _mm512_loadu_ps(&xre[i+112]);
                             zmm18  = _mm512_loadu_ps(&yre[i+112]);
                             zmm19  = _mm512_loadu_ps(&xim[i+112]);
                             zmm20  = _mm512_loadu_ps(&yim[i+112]);
                             zmm21  = _mm512_fmadd_ps(zmm18,zmm18,
                                                          _mm512_mul_ps(zmm20,zmm20)); // den
                             zmm22  = _mm512_fmsub_ps(zmm17,zmm18,
                                                          _mm512_mul_ps(zmm19,zmm20));// rep
                             redr[7] = _mm512_add_ps(redr[7],_mm512_div_ps(zmm22,zmm21)); 
                             zmm23   = _mm512_fmadd_ps(zmm19,zmm18,
                                                          _mm512_mul_ps(zmm17,zmm20)); // imp
                             redi[7] = _mm512_add_ps(redi[7],_mm512_div_ps(zmm23,zmm21)); 
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm24  = _mm512_loadu_ps(&xre[i+128]);
                             zmm25  = _mm512_loadu_ps(&yre[i+128]);
                             zmm26  = _mm512_loadu_ps(&xim[i+128]);
                             zmm27  = _mm512_loadu_ps(&yim[i+128]);
                             zmm28  = _mm512_fmadd_ps(zmm25,zmm25,
                                                          _mm512_mul_ps(zmm27,zmm27)); // den
                             zmm29  = _mm512_fmsub_ps(zmm24,zmm25,
                                                          _mm512_mul_ps(zmm26,zmm27));// rep
                             redr[8] = _mm512_add_ps(redr[8],_mm512_div_ps(zmm29,zmm28)); 
                             zmm30   = _mm512_fmadd_ps(zmm26,zmm25,
                                                          _mm512_mul_ps(zmm24,zmm27)); // imp
                             redi[8] = _mm512_add_ps(redi[8],_mm512_div_ps(zmm30,zmm28));
                             zmm30   = _mm512_loadu_ps(&xre[i+144]);
                             zmm31   = _mm512_loadu_ps(&yre[i+144]);
                             zmm0    = _mm512_loadu_ps(&xim[i+144]);
                             zmm1    = _mm512_loadu_ps(&yim[i+144]);
                             zmm2    = _mm512_fmadd_ps(zmm31,zmm31,
                                                          _mm512_mul_ps(zmm1,zmm1)); // den
                             zmm3   = _mm512_fmsub_ps(zmm30,zmm31,
                                                          _mm512_mul_ps(zmm0,zmm1));// rep
                             redr[9] = _mm512_add_ps(redr[9],_mm512_div_ps(zmm3,zmm2)); 
                             zmm4    = _mm512_fmadd_ps(zmm0,zmm31,
                                                          _mm512_mul_ps(zmm30,zmm1)); // imp
                             redi[9] = _mm512_add_ps(redi[9],_mm512_div_ps(zmm4,zmm2));  
                        }

                             redr[0] = _mm512_add_ps(redr[0],redr[5]);
                             redi[0] = _mm512_add_ps(redi[0],redi[5]);
                             redr[1] = _mm512_add_ps(redr[1],redr[6]);
                             redi[1] = _mm512_add_ps(redi[1],redi[6]);
                             redr[2] = _mm512_add_ps(redr[2],redr[7]);
                             redi[2] = _mm512_add_ps(redi[2],redi[7]);
                             redr[3] = _mm512_add_ps(redr[3],redr[8]);
                             redi[3] = _mm512_add_ps(redi[3],redi[8]);
                             redr[4] = _mm512_add_ps(redr[4],redr[9]);
                             redi[4] = _mm512_add_ps(redi[4],redi[9]); 

                         for(; (i+79) < n; i += 80) {
                             zmm0    = _mm512_loadu_ps(&xre[i+0]);
                             zmm1    = _mm512_loadu_ps(&yre[i+0]);
                             zmm2    = _mm512_loadu_ps(&xim[i+0]);
                             zmm3    = _mm512_loadu_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4));
                             zmm7    = _mm512_loadu_ps(&xre[i+16]);
                             zmm8    = _mm512_loadu_ps(&yre[i+16]);
                             zmm9    = _mm512_loadu_ps(&xim[i+16]);
                             zmm10   = _mm512_loadu_ps(&yim[i+16]);
                             zmm11   = _mm512_fmadd_ps(zmm8,zmm8,
                                                          _mm512_mul_ps(zmm10,zmm10)); // den
                             zmm12   = _mm512_fmsub_ps(zmm7,zmm8,
                                                          _mm512_mul_ps(zmm9,zmm10));// rep
                             redr[1] = _mm512_add_ps(redr[1],_mm512_div_ps(zmm12,zmm11)); 
                             zmm13   = _mm512_fmadd_ps(zmm9,zmm8,
                                                          _mm512_mul_ps(zmm7,zmm10)); // imp
                             redi[1] = _mm512_add_ps(redi[1],_mm512_div_ps(zmm13,zmm11));
                             zmm14   = _mm512_loadu_ps(&xre[i+32]);
                             zmm15   = _mm512_loadu_ps(&yre[i+32]);
                             zmm16   = _mm512_loadu_ps(&xim[i+32]);
                             zmm17   = _mm512_loadu_ps(&yim[i+32]);
                             zmm18   = _mm512_fmadd_ps(zmm15,zmm15,
                                                          _mm512_mul_ps(zmm17,zmm17)); // den
                             zmm19   = _mm512_fmsub_ps(zmm14,zmm15,
                                                          _mm512_mul_ps(zmm16,zmm17));// rep
                             redr[2] = _mm512_add_ps(redr[2],_mm512_div_ps(zmm19,zmm18)); 
                             zmm20   = _mm512_fmadd_ps(zmm16,zmm15,
                                                          _mm512_mul_ps(zmm14,zmm17)); // imp
                             redi[2] = _mm512_add_ps(redi[2],_mm512_div_ps(zmm20,zmm18));  
                             zmm21   = _mm512_loadu_ps(&xre[i+48]);
                             zmm22   = _mm512_loadu_ps(&yre[i+48]);
                             zmm23   = _mm512_loadu_ps(&xim[i+48]);
                             zmm24   = _mm512_loadu_ps(&yim[i+48]);
                             zmm25   = _mm512_fmadd_ps(zmm22,zmm22,
                                                          _mm512_mul_ps(zmm24,zmm24)); // den
                             zmm26   = _mm512_fmsub_ps(zmm21,zmm22,
                                                          _mm512_mul_ps(zmm23,zmm24));// rep
                             redr[3] = _mm512_add_ps(redr[3],_mm512_div_ps(zmm26,zmm25)); 
                             zmm27   = _mm512_fmadd_ps(zmm23,zmm22,
                                                          _mm512_mul_ps(zmm21,zmm24)); // imp
                             redi[3] = _mm512_add_ps(redi[3],_mm512_div_ps(zmm27,zmm25)); 
                             zmm28   = _mm512_loadu_ps(&xre[i+64]);
                             zmm29   = _mm512_loadu_ps(&yre[i+64]);
                             zmm30   = _mm512_loadu_ps(&xim[i+64]);
                             zmm31   = _mm512_loadu_ps(&yim[i+64]);
                             zmm0    = _mm512_fmadd_ps(zmm29,zmm29,
                                                          _mm512_mul_ps(zmm31,zmm31)); // den
                             zmm1    = _mm512_fmsub_ps(zmm28,zmm29,
                                                          _mm512_mul_ps(zmm30,zmm31));// rep
                             redr[4] = _mm512_add_ps(redr[4],_mm512_div_ps(zmm1,zmm0)); 
                             zmm2    = _mm512_fmadd_ps(zmm30,zmm29,
                                                          _mm512_mul_ps(zmm28,zmm31)); // imp
                             redi[4] = _mm512_add_ps(redi[4],_mm512_div_ps(zmm2,zmm0)); 
                        }

                          redr[0] = _mm512_add_ps(redr[0],redr[2]);
                          redi[0] = _mm512_add_ps(redi[0],redi[2]);
                          redr[1] = _mm512_add_ps(redr[1],redr[3]);
                          redi[1] = _mm512_add_ps(redi[1],redi[3]);
                          redr[0] = _mm512_add_ps(redr[0],redr[4]);
                          redi[0] = _mm512_add_ps(redi[0],redi[4]);

                           for(; (i+31) < n; i += 32) {
                             zmm0    = _mm512_loadu_ps(&xre[i+0]);
                             zmm1    = _mm512_loadu_ps(&yre[i+0]);
                             zmm2    = _mm512_loadu_ps(&xim[i+0]);
                             zmm3    = _mm512_loadu_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4));
                             zmm7    = _mm512_loadu_ps(&xre[i+16]);
                             zmm8    = _mm512_loadu_ps(&yre[i+16]);
                             zmm9    = _mm512_loadu_ps(&xim[i+16]);
                             zmm10   = _mm512_loadu_ps(&yim[i+16]);
                             zmm11   = _mm512_fmadd_ps(zmm8,zmm8,
                                                          _mm512_mul_ps(zmm10,zmm10)); // den
                             zmm12   = _mm512_fmsub_ps(zmm7,zmm8,
                                                          _mm512_mul_ps(zmm9,zmm10));// rep
                             redr[1] = _mm512_add_ps(redr[1],_mm512_div_ps(zmm12,zmm11)); 
                             zmm13   = _mm512_fmadd_ps(zmm9,zmm8,
                                                          _mm512_mul_ps(zmm7,zmm10)); // imp
                             redi[1] = _mm512_add_ps(redi[1],_mm512_div_ps(zmm13,zmm11));
                        }

                           redr[0] = _mm512_add_ps(redr[0],redr[1]);
                           redi[0] = _mm512_add_ps(redi[0],redi[1]);  

                         for(; (i+15) < n; i += 16) {
                             zmm0    = _mm512_loadu_ps(&xre[i+0]);
                             zmm1    = _mm512_loadu_ps(&yre[i+0]);
                             zmm2    = _mm512_loadu_ps(&xim[i+0]);
                             zmm3    = _mm512_loadu_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4)); 
                        }

                         for(; (i+1) < n; i += 1) {
                             const float xr = xre[i];
                             const float yr = yre[i];
                             const float xi = xim[i];
                             const float yi = yim[i];
                             const float den= (yr*yr)+(yi*yi);
                             re             +=((xr*yr)-(xi*yi))/den;
                             im             +=((xi*yr)+(xr*yi))/den;
                       }
                       
                       re   += _mm512_reduce_add_ps(redr[0]);
                       *mre =  re*invN;
                       im   += _mm512_reduce_add_ps(redi[0]);
                       *mim =  im*invN; 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void cmean_quot_u10x_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict __ATTR_ALIGN__(64) mre,
                                                  float * __restrict __ATTR_ALIGN__(64) mim,
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
                         __ATTR_ALIGN__(64) __m512 redr[10] = {_mm512_setzero_ps()};
                         __ATTR_ALIGN__(64) __m512 redi[10] = {_mm512_setzero_ps()};
                         register float re,im;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         re = 0.0f;
                         im = 0.0f;

                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0    = _mm512_load_ps(&xre[i+0]);
                             zmm1    = _mm512_load_ps(&yre[i+0]);
                             zmm2    = _mm512_load_ps(&xim[i+0]);
                             zmm3    = _mm512_load_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4));
                             zmm7    = _mm512_load_ps(&xre[i+16]);
                             zmm8    = _mm512_load_ps(&yre[i+16]);
                             zmm9    = _mm512_load_ps(&xim[i+16]);
                             zmm10   = _mm512_load_ps(&yim[i+16]);
                             zmm11   = _mm512_fmadd_ps(zmm8,zmm8,
                                                          _mm512_mul_ps(zmm10,zmm10)); // den
                             zmm12   = _mm512_fmsub_ps(zmm7,zmm8,
                                                          _mm512_mul_ps(zmm9,zmm10));// rep
                             redr[1] = _mm512_add_ps(redr[1],_mm512_div_ps(zmm12,zmm11)); 
                             zmm13   = _mm512_fmadd_ps(zmm9,zmm8,
                                                          _mm512_mul_ps(zmm7,zmm10)); // imp
                             redi[1] = _mm512_add_ps(redi[1],_mm512_div_ps(zmm13,zmm11));
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm14   = _mm512_load_ps(&xre[i+32]);
                             zmm15   = _mm512_load_ps(&yre[i+32]);
                             zmm16   = _mm512_load_ps(&xim[i+32]);
                             zmm17   = _mm512_load_ps(&yim[i+32]);
                             zmm18   = _mm512_fmadd_ps(zmm15,zmm15,
                                                          _mm512_mul_ps(zmm17,zmm17)); // den
                             zmm19   = _mm512_fmsub_ps(zmm14,zmm15,
                                                          _mm512_mul_ps(zmm16,zmm17));// rep
                             redr[2] = _mm512_add_ps(redr[2],_mm512_div_ps(zmm19,zmm18)); 
                             zmm20   = _mm512_fmadd_ps(zmm16,zmm15,
                                                          _mm512_mul_ps(zmm14,zmm17)); // imp
                             redi[2] = _mm512_add_ps(redi[2],_mm512_div_ps(zmm20,zmm18));  
                             zmm21   = _mm512_load_ps(&xre[i+48]);
                             zmm22   = _mm512_load_ps(&yre[i+48]);
                             zmm23   = _mm512_load_ps(&xim[i+48]);
                             zmm24   = _mm512_load_ps(&yim[i+48]);
                             zmm25   = _mm512_fmadd_ps(zmm22,zmm22,
                                                          _mm512_mul_ps(zmm24,zmm24)); // den
                             zmm26   = _mm512_fmsub_ps(zmm21,zmm22,
                                                          _mm512_mul_ps(zmm23,zmm24));// rep
                             redr[3] = _mm512_add_ps(redr[3],_mm512_div_ps(zmm26,zmm25)); 
                             zmm27   = _mm512_fmadd_ps(zmm23,zmm22,
                                                          _mm512_mul_ps(zmm21,zmm24)); // imp
                             redi[3] = _mm512_add_ps(redi[3],_mm512_div_ps(zmm27,zmm25)); 
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0); 
                             zmm28   = _mm512_load_ps(&xre[i+64]);
                             zmm29   = _mm512_load_ps(&yre[i+64]);
                             zmm30   = _mm512_load_ps(&xim[i+64]);
                             zmm31   = _mm512_load_ps(&yim[i+64]);
                             zmm0    = _mm512_fmadd_ps(zmm29,zmm29,
                                                          _mm512_mul_ps(zmm31,zmm31)); // den
                             zmm1    = _mm512_fmsub_ps(zmm28,zmm29,
                                                          _mm512_mul_ps(zmm30,zmm31));// rep
                             redr[4] = _mm512_add_ps(redr[4],_mm512_div_ps(zmm1,zmm0)); 
                             zmm2    = _mm512_fmadd_ps(zmm30,zmm29,
                                                          _mm512_mul_ps(zmm28,zmm31)); // imp
                             redi[4] = _mm512_add_ps(redi[4],_mm512_div_ps(zmm2,zmm0));   
                             zmm3   = _mm512_load_ps(&xre[i+80]);
                             zmm4   = _mm512_load_ps(&yre[i+80]);
                             zmm5   = _mm512_load_ps(&xim[i+80]);
                             zmm6   = _mm512_load_ps(&yim[i+80]);
                             zmm7   = _mm512_fmadd_ps(zmm4,zmm4,
                                                          _mm512_mul_ps(zmm6,zmm6)); // den
                             zmm8   = _mm512_fmsub_ps(zmm3,zmm4,
                                                          _mm512_mul_ps(zmm5,zmm6));// rep
                             redr[5]= _mm512_add_ps(redr[5],_mm512_div_ps(zmm8,zmm7)); 
                             zmm9   = _mm512_fmadd_ps(zmm5,zmm4,
                                                          _mm512_mul_ps(zmm3,zmm6)); // imp
                             redi[5]= _mm512_add_ps(redi[5],_mm512_div_ps(zmm9,zmm7)); 
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm10  = _mm512_load_ps(&xre[i+96]);
                             zmm11  = _mm512_load_ps(&yre[i+96]);
                             zmm12  = _mm512_load_ps(&xim[i+96]);
                             zmm13  = _mm512_load_ps(&yim[i+96]);
                             zmm14  = _mm512_fmadd_ps(zmm11,zmm11,
                                                          _mm512_mul_ps(zmm13,zmm13)); // den
                             zmm15   = _mm512_fmsub_ps(zmm10,zmm11,
                                                          _mm512_mul_ps(zmm12,zmm13));// rep
                             redr[6] = _mm512_add_ps(redr[6],_mm512_div_ps(zmm15,zmm14)); 
                             zmm16    = _mm512_fmadd_ps(zmm12,zmm11,
                                                          _mm512_mul_ps(zmm10,zmm13)); // imp
                             redi[6] = _mm512_add_ps(redi[6],_mm512_div_ps(zmm16,zmm14)); 
                             zmm17  = _mm512_load_ps(&xre[i+112]);
                             zmm18  = _mm512_load_ps(&yre[i+112]);
                             zmm19  = _mm512_load_ps(&xim[i+112]);
                             zmm20  = _mm512_load_ps(&yim[i+112]);
                             zmm21  = _mm512_fmadd_ps(zmm18,zmm18,
                                                          _mm512_mul_ps(zmm20,zmm20)); // den
                             zmm22  = _mm512_fmsub_ps(zmm17,zmm18,
                                                          _mm512_mul_ps(zmm19,zmm20));// rep
                             redr[7] = _mm512_add_ps(redr[7],_mm512_div_ps(zmm22,zmm21)); 
                             zmm23   = _mm512_fmadd_ps(zmm19,zmm18,
                                                          _mm512_mul_ps(zmm17,zmm20)); // imp
                             redi[7] = _mm512_add_ps(redi[7],_mm512_div_ps(zmm23,zmm21)); 
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm24  = _mm512_load_ps(&xre[i+128]);
                             zmm25  = _mm512_load_ps(&yre[i+128]);
                             zmm26  = _mm512_load_ps(&xim[i+128]);
                             zmm27  = _mm512_load_ps(&yim[i+128]);
                             zmm28  = _mm512_fmadd_ps(zmm25,zmm25,
                                                          _mm512_mul_ps(zmm27,zmm27)); // den
                             zmm29  = _mm512_fmsub_ps(zmm24,zmm25,
                                                          _mm512_mul_ps(zmm26,zmm27));// rep
                             redr[8] = _mm512_add_ps(redr[8],_mm512_div_ps(zmm29,zmm28)); 
                             zmm30   = _mm512_fmadd_ps(zmm26,zmm25,
                                                          _mm512_mul_ps(zmm24,zmm27)); // imp
                             redi[8] = _mm512_add_ps(redi[8],_mm512_div_ps(zmm30,zmm28));
                             zmm30   = _mm512_load_ps(&xre[i+144]);
                             zmm31   = _mm512_load_ps(&yre[i+144]);
                             zmm0    = _mm512_load_ps(&xim[i+144]);
                             zmm1    = _mm512_load_ps(&yim[i+144]);
                             zmm2    = _mm512_fmadd_ps(zmm31,zmm31,
                                                          _mm512_mul_ps(zmm1,zmm1)); // den
                             zmm3   = _mm512_fmsub_ps(zmm30,zmm31,
                                                          _mm512_mul_ps(zmm0,zmm1));// rep
                             redr[9] = _mm512_add_ps(redr[9],_mm512_div_ps(zmm3,zmm2)); 
                             zmm4    = _mm512_fmadd_ps(zmm0,zmm31,
                                                          _mm512_mul_ps(zmm30,zmm1)); // imp
                             redi[9] = _mm512_add_ps(redi[9],_mm512_div_ps(zmm4,zmm2));  
                        }

                             redr[0] = _mm512_add_ps(redr[0],redr[5]);
                             redi[0] = _mm512_add_ps(redi[0],redi[5]);
                             redr[1] = _mm512_add_ps(redr[1],redr[6]);
                             redi[1] = _mm512_add_ps(redi[1],redi[6]);
                             redr[2] = _mm512_add_ps(redr[2],redr[7]);
                             redi[2] = _mm512_add_ps(redi[2],redi[7]);
                             redr[3] = _mm512_add_ps(redr[3],redr[8]);
                             redi[3] = _mm512_add_ps(redi[3],redi[8]);
                             redr[4] = _mm512_add_ps(redr[4],redr[9]);
                             redi[4] = _mm512_add_ps(redi[4],redi[9]); 

                         for(; (i+79) < n; i += 80) {
                             zmm0    = _mm512_load_ps(&xre[i+0]);
                             zmm1    = _mm512_load_ps(&yre[i+0]);
                             zmm2    = _mm512_load_ps(&xim[i+0]);
                             zmm3    = _mm512_load_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4));
                             zmm7    = _mm512_load_ps(&xre[i+16]);
                             zmm8    = _mm512_load_ps(&yre[i+16]);
                             zmm9    = _mm512_load_ps(&xim[i+16]);
                             zmm10   = _mm512_load_ps(&yim[i+16]);
                             zmm11   = _mm512_fmadd_ps(zmm8,zmm8,
                                                          _mm512_mul_ps(zmm10,zmm10)); // den
                             zmm12   = _mm512_fmsub_ps(zmm7,zmm8,
                                                          _mm512_mul_ps(zmm9,zmm10));// rep
                             redr[1] = _mm512_add_ps(redr[1],_mm512_div_ps(zmm12,zmm11)); 
                             zmm13   = _mm512_fmadd_ps(zmm9,zmm8,
                                                          _mm512_mul_ps(zmm7,zmm10)); // imp
                             redi[1] = _mm512_add_ps(redi[1],_mm512_div_ps(zmm13,zmm11));
                             zmm14   = _mm512_load_ps(&xre[i+32]);
                             zmm15   = _mm512_load_ps(&yre[i+32]);
                             zmm16   = _mm512_load_ps(&xim[i+32]);
                             zmm17   = _mm512_load_ps(&yim[i+32]);
                             zmm18   = _mm512_fmadd_ps(zmm15,zmm15,
                                                          _mm512_mul_ps(zmm17,zmm17)); // den
                             zmm19   = _mm512_fmsub_ps(zmm14,zmm15,
                                                          _mm512_mul_ps(zmm16,zmm17));// rep
                             redr[2] = _mm512_add_ps(redr[2],_mm512_div_ps(zmm19,zmm18)); 
                             zmm20   = _mm512_fmadd_ps(zmm16,zmm15,
                                                          _mm512_mul_ps(zmm14,zmm17)); // imp
                             redi[2] = _mm512_add_ps(redi[2],_mm512_div_ps(zmm20,zmm18));  
                             zmm21   = _mm512_load_ps(&xre[i+48]);
                             zmm22   = _mm512_load_ps(&yre[i+48]);
                             zmm23   = _mm512_load_ps(&xim[i+48]);
                             zmm24   = _mm512_load_ps(&yim[i+48]);
                             zmm25   = _mm512_fmadd_ps(zmm22,zmm22,
                                                          _mm512_mul_ps(zmm24,zmm24)); // den
                             zmm26   = _mm512_fmsub_ps(zmm21,zmm22,
                                                          _mm512_mul_ps(zmm23,zmm24));// rep
                             redr[3] = _mm512_add_ps(redr[3],_mm512_div_ps(zmm26,zmm25)); 
                             zmm27   = _mm512_fmadd_ps(zmm23,zmm22,
                                                          _mm512_mul_ps(zmm21,zmm24)); // imp
                             redi[3] = _mm512_add_ps(redi[3],_mm512_div_ps(zmm27,zmm25)); 
                             zmm28   = _mm512_load_ps(&xre[i+64]);
                             zmm29   = _mm512_load_ps(&yre[i+64]);
                             zmm30   = _mm512_load_ps(&xim[i+64]);
                             zmm31   = _mm512_load_ps(&yim[i+64]);
                             zmm0    = _mm512_fmadd_ps(zmm29,zmm29,
                                                          _mm512_mul_ps(zmm31,zmm31)); // den
                             zmm1    = _mm512_fmsub_ps(zmm28,zmm29,
                                                          _mm512_mul_ps(zmm30,zmm31));// rep
                             redr[4] = _mm512_add_ps(redr[4],_mm512_div_ps(zmm1,zmm0)); 
                             zmm2    = _mm512_fmadd_ps(zmm30,zmm29,
                                                          _mm512_mul_ps(zmm28,zmm31)); // imp
                             redi[4] = _mm512_add_ps(redi[4],_mm512_div_ps(zmm2,zmm0)); 
                        }

                          redr[0] = _mm512_add_ps(redr[0],redr[2]);
                          redi[0] = _mm512_add_ps(redi[0],redi[2]);
                          redr[1] = _mm512_add_ps(redr[1],redr[3]);
                          redi[1] = _mm512_add_ps(redi[1],redi[3]);
                          redr[0] = _mm512_add_ps(redr[0],redr[4]);
                          redi[0] = _mm512_add_ps(redi[0],redi[4]);

                           for(; (i+31) < n; i += 32) {
                             zmm0    = _mm512_load_ps(&xre[i+0]);
                             zmm1    = _mm512_load_ps(&yre[i+0]);
                             zmm2    = _mm512_load_ps(&xim[i+0]);
                             zmm3    = _mm512_load_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4));
                             zmm7    = _mm512_load_ps(&xre[i+16]);
                             zmm8    = _mm512_load_ps(&yre[i+16]);
                             zmm9    = _mm512_load_ps(&xim[i+16]);
                             zmm10   = _mm512_load_ps(&yim[i+16]);
                             zmm11   = _mm512_fmadd_ps(zmm8,zmm8,
                                                          _mm512_mul_ps(zmm10,zmm10)); // den
                             zmm12   = _mm512_fmsub_ps(zmm7,zmm8,
                                                          _mm512_mul_ps(zmm9,zmm10));// rep
                             redr[1] = _mm512_add_ps(redr[1],_mm512_div_ps(zmm12,zmm11)); 
                             zmm13   = _mm512_fmadd_ps(zmm9,zmm8,
                                                          _mm512_mul_ps(zmm7,zmm10)); // imp
                             redi[1] = _mm512_add_ps(redi[1],_mm512_div_ps(zmm13,zmm11));
                        }

                           redr[0] = _mm512_add_ps(redr[0],redr[1]);
                           redi[0] = _mm512_add_ps(redi[0],redi[1]);  

                         for(; (i+15) < n; i += 16) {
                             zmm0    = _mm512_load_ps(&xre[i+0]);
                             zmm1    = _mm512_load_ps(&yre[i+0]);
                             zmm2    = _mm512_load_ps(&xim[i+0]);
                             zmm3    = _mm512_load_ps(&yim[i+0]);
                             zmm4    = _mm512_fmadd_ps(zmm1,zmm1,
                                                          _mm512_mul_ps(zmm3,zmm3)); // den
                             zmm5    = _mm512_fmsub_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3));// rep
                             redr[0] = _mm512_add_ps(redr[0],_mm512_div_ps(zmm5,zmm4)); 
                             zmm6    = _mm512_fmadd_ps(zmm2,zmm1,
                                                          _mm512_mul_ps(zmm0,zmm3)); // imp
                             redi[0] = _mm512_add_ps(redi[0],_mm512_div_ps(zmm6,zmm4)); 
                        }

                         for(; (i+1) < n; i += 1) {
                             const float xr = xre[i];
                             const float yr = yre[i];
                             const float xi = xim[i];
                             const float yi = yim[i];
                             const float den= (yr*yr)+(yi*yi);
                             re             +=((xr*yr)-(xi*yi))/den;
                             im             +=((xi*yr)+(xr*yi))/den;
                       }
                       
                       re   += _mm512_reduce_add_ps(redr[0]);
                       *mre =  re*invN;
                       im   += _mm512_reduce_add_ps(redi[0]);
                       *mim =  im*invN; 
               }




     }// math



} //gms





































#endif /*__GMS_CMEAN_QUOT_ZMM16R4_HPP__*/
