

#ifndef __GMS_CARITHM_MEAN_ZMM16R4_HPP__
#define __GMS_CARITHM_MEAN_ZMM16R4_HPP__ 130420230909


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

    const unsigned int GMS_CARITHM_MEAN_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CARITHM_MEAN_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CARITHM_MEAN_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CARITHM_MEAN_ZMM16R4_FULLVER =
      1000U*GMS_CARITHM_MEAN_ZMM16R4_MAJOR+
      100U*GMS_CARITHM_MEAN_ZMM16R4_MINOR+
      10U*GMS_CARITHM_MEAN_ZMM16R4_MICRO;
    const char * const GMS_CARITHM_MEAN_ZMM16R4_CREATION_DATE = "13-04-2023 09:09 AM +00200 ( THR 13 APR 2023 GMT+2)";
    const char * const GMS_CARITHM_MEAN_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CARITHM_MEAN_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CARITHM_MEAN_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex arithmetic mean operations."

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
                   void cmean_arithm_u10x_zmm16r4_u(const float * __restrict xre,
                                                    const float * __restrict xim,
                                                    float * __restrict mre,
                                                    float * __restrict mim,
                                                    const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                         register __m512 zmm0,zmm1,zmm2,zmm3;
                         register __m512 zmm4,zmm5,zmm6,zmm7;
                         register __m512 zmm8,zmm9,zmm10,zmm11;
                         register __m512 zmm12,zmm13,zmm14,zmm15;
                         register __m512 zmm16,zmm17,zmm18,zmm19;
                         __m512 redr[10] = {_mm512_setzero_ps()};
                         __m512 redi[10] = {_mm512_setzero_ps()};
                         register float re,im;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         re = 0.0f;
                         im = 0.0f;  
                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                            zmm0    = _mm512_loadu_ps(&xre[i+0]);
                            zmm1    = _mm512_loadu_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                            zmm2    = _mm512_loadu_ps(&xre[i+16]);
                            zmm3    = _mm512_loadu_ps(&xim[i+16]);
                            redr[1] = _mm512_add_pd(redr[1],zmm2);
                            redi[1] = _mm512_add_ps(redi[1],zmm3);
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                            zmm4    = _mm512_loadu_ps(&xre[i+32]);
                            zmm5    = _mm512_loadu_ps(&xim[i+32]);
                            redr[2] = _mm512_add_pd(redr[2],zmm4);
                            redi[2] = _mm512_add_ps(redi[2],zmm5);
                            zmm6    = _mm512_loadu_ps(&xre[i+48]);
                            zmm7    = _mm512_loadu_ps(&xim[i+48]);
                            redr[3] = _mm512_add_pd(redr[3],zmm6);
                            redi[3] = _mm512_add_ps(redi[3],zmm7);
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                            zmm8    = _mm512_loadu_ps(&xre[i+64]);
                            zmm9    = _mm512_loadu_ps(&xim[i+64]);
                            redr[4] = _mm512_add_pd(redr[4],zmm8);
                            redi[4] = _mm512_add_ps(redi[4],zmm9);
                            zmm10   = _mm512_loadu_ps(&xre[i+80]);
                            zmm11   = _mm512_loadu_ps(&xim[i+80]);
                            redr[5] = _mm512_add_pd(redr[5],zmm10);
                            redi[5] = _mm512_add_ps(redi[5],zmm11);
                            _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                            zmm12   = _mm512_loadu_ps(&xre[i+96]);
                            zmm13   = _mm512_loadu_ps(&xim[i+96]);
                            redr[6] = _mm512_add_pd(redr[6],zmm12);
                            redi[6] = _mm512_add_ps(redi[6],zmm13);
                            zmm14   = _mm512_loadu_ps(&xre[i+112]);
                            zmm15   = _mm512_loadu_ps(&xim[i+112]);
                            redr[7] = _mm512_add_pd(redr[7],zmm14);
                            redi[7] = _mm512_add_ps(redi[7],zmm15);
                            _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0); 
                            zmm16   = _mm512_loadu_ps(&xre[i+128]);
                            zmm17   = _mm512_loadu_ps(&xim[i+128]);
                            redr[8] = _mm512_add_pd(redr[8],zmm16);
                            redi[8] = _mm512_add_ps(redi[8],zmm17);
                            zmm18   = _mm512_loadu_ps(&xre[i+144]);
                            zmm19   = _mm512_loadu_ps(&xim[i+144]);
                            redr[9] = _mm512_add_pd(redr[9],zmm18);
                            redi[9] = _mm512_add_ps(redi[9],zmm19);
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
                            zmm1    = _mm512_loadu_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                            zmm2    = _mm512_loadu_ps(&xre[i+16]);
                            zmm3    = _mm512_loadu_ps(&xim[i+16]);
                            redr[1] = _mm512_add_pd(redr[1],zmm2);
                            redi[1] = _mm512_add_ps(redi[1],zmm3);
                            zmm4    = _mm512_loadu_ps(&xre[i+32]);
                            zmm5    = _mm512_loadu_ps(&xim[i+32]);
                            redr[2] = _mm512_add_pd(redr[2],zmm4);
                            redi[2] = _mm512_add_ps(redi[2],zmm5);
                            zmm6    = _mm512_loadu_ps(&xre[i+48]);
                            zmm7    = _mm512_loadu_ps(&xim[i+48]);
                            redr[3] = _mm512_add_pd(redr[3],zmm6);
                            redi[3] = _mm512_add_ps(redi[3],zmm7);
                            zmm8    = _mm512_loadu_ps(&xre[i+64]);
                            zmm9    = _mm512_loadu_ps(&xim[i+64]);
                            redr[4] = _mm512_add_pd(redr[4],zmm8);
                            redi[4] = _mm512_add_ps(redi[4],zmm9); 
                      }

                             redr[0] = _mm512_add_ps(redr[0],redr[2]);
                             redi[0] = _mm512_add_ps(redi[0],redi[2]);
                             redr[1] = _mm512_add_ps(redr[1],redr[3]);
                             redi[1] = _mm512_add_ps(redi[1],redi[3]);
                             redr[0] = _mm512_add_ps(redr[0],redr[4]);
                             redi[0] = _mm512_add_ps(redi[0],redi[4]);

                       for(; (i+31) < n; i += 32) {
                            zmm0    = _mm512_loadu_ps(&xre[i+0]);
                            zmm1    = _mm512_loadu_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                            zmm2    = _mm512_loadu_ps(&xre[i+16]);
                            zmm3    = _mm512_loadu_ps(&xim[i+16]);
                            redr[1] = _mm512_add_pd(redr[1],zmm2);
                            redi[1] = _mm512_add_ps(redi[1],zmm3); 
                      }

                              redr[0] = _mm512_add_ps(redr[0],redr[1]);
                              redi[0] = _mm512_add_ps(redi[0],redi[1]); 

                      for(; (i+15) < n; i += 16) {
                            zmm0    = _mm512_loadu_ps(&xre[i+0]);
                            zmm1    = _mm512_loadu_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                     } 

                      for(; (i+0) < n; i += 1) {
                            const float xr = xre[i];
                            const float xi = xim[i];
                            re += xr;
                            im += xi;
                     }

                      re += _mm512_reduce_add_ps(redr[0]);
                      *mre = re*invN;
                      im += _mm512_reduce_add_ps(redi[0]);
                      *mim = im*invN;
              }   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline     
                   void cmean_arithm_u10x_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                    const float * __restrict __ATTR_ALIGN__(64) xim,
                                                    float * __restrict mre,
                                                    float * __restrict mim,
                                                    const int32_t n) {

                        if(__builtin_expect(0==n,0)) {return;}
                         register __m512 zmm0,zmm1,zmm2,zmm3;
                         register __m512 zmm4,zmm5,zmm6,zmm7;
                         register __m512 zmm8,zmm9,zmm10,zmm11;
                         register __m512 zmm12,zmm13,zmm14,zmm15;
                         register __m512 zmm16,zmm17,zmm18,zmm19;
                         __ATTR_ALIGN__(64)  __m512 redr[10] = {_mm512_setzero_ps()};
                         __ATTR_ALIGN__(64)  __m512 redi[10] = {_mm512_setzero_ps()};
                         register float re,im;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         re = 0.0f;
                         im = 0.0f;  
                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                            zmm0    = _mm512_load_ps(&xre[i+0]);
                            zmm1    = _mm512_load_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                            zmm2    = _mm512_load_ps(&xre[i+16]);
                            zmm3    = _mm512_load_ps(&xim[i+16]);
                            redr[1] = _mm512_add_pd(redr[1],zmm2);
                            redi[1] = _mm512_add_ps(redi[1],zmm3);
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                            zmm4    = _mm512_load_ps(&xre[i+32]);
                            zmm5    = _mm512_load_ps(&xim[i+32]);
                            redr[2] = _mm512_add_pd(redr[2],zmm4);
                            redi[2] = _mm512_add_ps(redi[2],zmm5);
                            zmm6    = _mm512_load_ps(&xre[i+48]);
                            zmm7    = _mm512_load_ps(&xim[i+48]);
                            redr[3] = _mm512_add_pd(redr[3],zmm6);
                            redi[3] = _mm512_add_ps(redi[3],zmm7);
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                            zmm8    = _mm512_load_ps(&xre[i+64]);
                            zmm9    = _mm512_load_ps(&xim[i+64]);
                            redr[4] = _mm512_add_pd(redr[4],zmm8);
                            redi[4] = _mm512_add_ps(redi[4],zmm9);
                            zmm10   = _mm512_load_ps(&xre[i+80]);
                            zmm11   = _mm512_load_ps(&xim[i+80]);
                            redr[5] = _mm512_add_pd(redr[5],zmm10);
                            redi[5] = _mm512_add_ps(redi[5],zmm11);
                            _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                            zmm12   = _mm512_load_ps(&xre[i+96]);
                            zmm13   = _mm512_load_ps(&xim[i+96]);
                            redr[6] = _mm512_add_pd(redr[6],zmm12);
                            redi[6] = _mm512_add_ps(redi[6],zmm13);
                            zmm14   = _mm512_load_ps(&xre[i+112]);
                            zmm15   = _mm512_load_ps(&xim[i+112]);
                            redr[7] = _mm512_add_pd(redr[7],zmm14);
                            redi[7] = _mm512_add_ps(redi[7],zmm15);
                            _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0); 
                            zmm16   = _mm512_load_ps(&xre[i+128]);
                            zmm17   = _mm512_load_ps(&xim[i+128]);
                            redr[8] = _mm512_add_pd(redr[8],zmm16);
                            redi[8] = _mm512_add_ps(redi[8],zmm17);
                            zmm18   = _mm512_load_ps(&xre[i+144]);
                            zmm19   = _mm512_load_ps(&xim[i+144]);
                            redr[9] = _mm512_add_pd(redr[9],zmm18);
                            redi[9] = _mm512_add_ps(redi[9],zmm19);
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
                            zmm1    = _mm512_load_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                            zmm2    = _mm512_load_ps(&xre[i+16]);
                            zmm3    = _mm512_load_ps(&xim[i+16]);
                            redr[1] = _mm512_add_pd(redr[1],zmm2);
                            redi[1] = _mm512_add_ps(redi[1],zmm3);
                            zmm4    = _mm512_load_ps(&xre[i+32]);
                            zmm5    = _mm512_load_ps(&xim[i+32]);
                            redr[2] = _mm512_add_pd(redr[2],zmm4);
                            redi[2] = _mm512_add_ps(redi[2],zmm5);
                            zmm6    = _mm512_load_ps(&xre[i+48]);
                            zmm7    = _mm512_load_ps(&xim[i+48]);
                            redr[3] = _mm512_add_pd(redr[3],zmm6);
                            redi[3] = _mm512_add_ps(redi[3],zmm7);
                            zmm8    = _mm512_load_ps(&xre[i+64]);
                            zmm9    = _mm512_load_ps(&xim[i+64]);
                            redr[4] = _mm512_add_pd(redr[4],zmm8);
                            redi[4] = _mm512_add_ps(redi[4],zmm9); 
                      }

                             redr[0] = _mm512_add_ps(redr[0],redr[2]);
                             redi[0] = _mm512_add_ps(redi[0],redi[2]);
                             redr[1] = _mm512_add_ps(redr[1],redr[3]);
                             redi[1] = _mm512_add_ps(redi[1],redi[3]);
                             redr[0] = _mm512_add_ps(redr[0],redr[4]);
                             redi[0] = _mm512_add_ps(redi[0],redi[4]);

                       for(; (i+31) < n; i += 32) {
                            zmm0    = _mm512_load_ps(&xre[i+0]);
                            zmm1    = _mm512_load_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                            zmm2    = _mm512_load_ps(&xre[i+16]);
                            zmm3    = _mm512_load_ps(&xim[i+16]);
                            redr[1] = _mm512_add_pd(redr[1],zmm2);
                            redi[1] = _mm512_add_ps(redi[1],zmm3); 
                      }

                              redr[0] = _mm512_add_ps(redr[0],redr[1]);
                              redi[0] = _mm512_add_ps(redi[0],redi[1]); 

                      for(; (i+15) < n; i += 16) {
                            zmm0    = _mm512_load_ps(&xre[i+0]);
                            zmm1    = _mm512_load_ps(&xim[i+0]);
                            redr[0] = _mm512_add_pd(redr[0],zmm0);
                            redi[0] = _mm512_add_ps(redi[0],zmm1);
                     } 

                      for(; (i+0) < n; i += 1) {
                            const float xr = xre[i];
                            const float xi = xim[i];
                            re += xr;
                            im += xi;
                     }

                      re += _mm512_reduce_add_ps(redr[0]);
                      *mre = re*invN;
                      im += _mm512_reduce_add_ps(redi[0]);
                      *mim = im*invN;
              }     
  


      } // math

 




} // gms
































#endif /*__GMS_CARITHM_MEAN_ZMM16R4_HPP__*/
