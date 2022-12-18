
#ifndef __GMS_CPOLAR_ZMM16R4_HPP__
#define __GMS_CPOLAR_ZMM16R4_HPP__ 181220221013


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

    const unsigned int GMS_CPOLAR_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CPOLAR_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CPOLAR_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CPOLAR_ZMM16R4_FULLVER =
      1000U*GMS_CPOLAR_ZMM16R4_MAJOR+
      100U*GMS_CPOLAR_ZMM16R4_MINOR+
      10U*GMS_CPOLAR_ZMM16R4_MICRO;
    const char * const GMS_CPOLAR_ZMM16R4_CREATION_DATE = "18-12-2022 10:13 AM +00200 (SUN 18 DEC 2022 GMT+2)";
    const char * const GMS_CPOLAR_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CPOLAR_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CPOLAR_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex polar form."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"



namespace  gms {

         namespace  math {

                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cpolarv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) rho,
                                                     const float * __restrict __ATTR_ALIGN__(64) tht,
                                                     float * __restrict __ATTR_ALIGN__(64) re,
                                                     float * __restrict __ATTR_ALIGN__(64) im,
                                                     const int32_t n) {

                          if(__builtin_expect(0==n,0)) { return;}
                          register vfloat zmm0,zmm1,zmm2,zmm3;
                          register vfloat zmm4,zmm5,zmm6,zmm7;
                          register vfloat zmm8,zmm9,zmm10,zmm11;
                          register vfloat zmm12,zmm13,zmm14,zmm15;
                          register vfloat zmm16,zmm17,zmm18,zmm19;
                          register vfloat zmm20,zmm21,zmm22,zmm23;
                          register vfloat zmm24,zmm25,zmm26,zmm27;
                          register vfloat zmm28,zmm29,zmm30,zmm31;
                          int32_t i;
                          
                          for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char*)&rho[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&tht[i+32],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&rho[i+0]);
                              zmm1 = _mm512_load_ps(&tht[i+0]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                              _mm512_store_ps(&re[i+0],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                              _mm512_store_ps(&im[i+0],zmm3);
                              zmm4 = _mm512_load_ps(&rho[i+16]);
                              zmm5 = _mm512_load_ps(&tht[i+16]);
                              zmm6 = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&re[i+16],zmm6);
                              zmm7 = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&im[i+16],zmm7);
                              _mm_prefetch((const char*)&rho[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&tht[i+64],_MM_HINT_T0);
                              zmm8 = _mm512_load_ps(&rho[i+32]);
                              zmm9 = _mm512_load_ps(&tht[i+32]);
                              zmm10= _mm512_mul_ps(zmm8,xcosf(zmm9));
                              _mm512_store_ps(&re[i+32],zmm10);
                              zmm11= _mm512_mul_ps(zmm8,xsinf(zmm9));
                              _mm512_store_ps(&im[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&rho[i+48]);
                              zmm13 = _mm512_load_ps(&tht[i+48]);
                              zmm14 = _mm512_mul_ps(zmm12,xcosf(zmm13));
                              _mm512_store_ps(&re[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(zmm12,xsinf(zmm13));
                              _mm512_store_ps(&im[i+48],zmm15);
                              _mm_prefetch((const char*)&rho[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&tht[i+96],_MM_HINT_T0);
                              zmm16 = _mm512_load_ps(&rho[i+64]);
                              zmm17 = _mm512_load_ps(&tht[i+64]);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&re[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&im[i+64],zmm19);
                              zmm20 = _mm512_load_ps(&rho[i+80]);
                              zmm21 = _mm512_load_ps(&tht[i+80]);
                              zmm22 = _mm512_mul_ps(zmm20,xcosf(zmm21));
                              _mm512_store_ps(&re[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(zmm20,xsinf(zmm21));
                              _mm512_store_ps(&im[i+80],zmm23);
                              _mm_prefetch((const char*)&rho[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&tht[i+128],_MM_HINT_T0);
                              zmm24 = _mm512_load_ps(&rho[i+96]);
                              zmm25 = _mm512_load_ps(&tht[i+96]);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_store_ps(&re[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_store_ps(&im[i+96],zmm27);
                              zmm28 = _mm512_load_ps(&rho[i+112]);
                              zmm29 = _mm512_load_ps(&tht[i+112]);
                              zmm30 = _mm512_mul_ps(zmm28,xcosf(zmm29));
                              _mm512_store_ps(&re[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(zmm28,xsinf(zmm29));
                              _mm512_store_ps(&im[i+112],zmm31);
                              _mm_prefetch((const char*)&rho[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&tht[i+160],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&rho[i+128]);
                              zmm1 = _mm512_load_ps(&tht[i+128]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1));
                              _mm512_store_ps(&re[i+128],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1));
                              _mm512_store_ps(&im[i+128],zmm3);
                              zmm4 = _mm512_load_ps(&rho[i+144]);
                              zmm5 = _mm512_load_ps(&tht[i+144]);
                              zmm6 = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&re[i+144],zmm6);
                              zmm7 = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&im[i+144],zmm7);
                         }

                         for(; (i+127) < n; i += 128) {
                               zmm0 = _mm512_load_ps(&rho[i+0]);
                              zmm1 = _mm512_load_ps(&tht[i+0]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                              _mm512_store_ps(&re[i+0],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                              _mm512_store_ps(&im[i+0],zmm3);
                              zmm4 = _mm512_load_ps(&rho[i+16]);
                              zmm5 = _mm512_load_ps(&tht[i+16]);
                              zmm6 = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&re[i+16],zmm6);
                              zmm7 = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&im[i+16],zmm7);
                              zmm8 = _mm512_load_ps(&rho[i+32]);
                              zmm9 = _mm512_load_ps(&tht[i+32]);
                              zmm10= _mm512_mul_ps(zmm8,xcosf(zmm9));
                              _mm512_store_ps(&re[i+32],zmm10);
                              zmm11= _mm512_mul_ps(zmm8,xsinf(zmm9));
                              _mm512_store_ps(&im[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&rho[i+48]);
                              zmm13 = _mm512_load_ps(&tht[i+48]);
                              zmm14 = _mm512_mul_ps(zmm12,xcosf(zmm13));
                              _mm512_store_ps(&re[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(zmm12,xsinf(zmm13));
                              _mm512_store_ps(&im[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&rho[i+64]);
                              zmm17 = _mm512_load_ps(&tht[i+64]);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&re[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&im[i+64],zmm19);
                              zmm20 = _mm512_load_ps(&rho[i+80]);
                              zmm21 = _mm512_load_ps(&tht[i+80]);
                              zmm22 = _mm512_mul_ps(zmm20,xcosf(zmm21));
                              _mm512_store_ps(&re[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(zmm20,xsinf(zmm21));
                              _mm512_store_ps(&im[i+80],zmm23);
                              zmm24 = _mm512_load_ps(&rho[i+96]);
                              zmm25 = _mm512_load_ps(&tht[i+96]);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_store_ps(&re[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_store_ps(&im[i+96],zmm27);
                              zmm28 = _mm512_load_ps(&rho[i+112]);
                              zmm29 = _mm512_load_ps(&tht[i+112]);
                              zmm30 = _mm512_mul_ps(zmm28,xcosf(zmm29));
                              _mm512_store_ps(&re[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(zmm28,xsinf(zmm29));
                              _mm512_store_ps(&im[i+112],zmm31);

                         }

                          for(; (i+79) < n; i += 80) {
                              zmm0 = _mm512_load_ps(&rho[i+0]);
                              zmm1 = _mm512_load_ps(&tht[i+0]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                              _mm512_store_ps(&re[i+0],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                              _mm512_store_ps(&im[i+0],zmm3);
                              zmm4 = _mm512_load_ps(&rho[i+16]);
                              zmm5 = _mm512_load_ps(&tht[i+16]);
                              zmm6 = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&re[i+16],zmm6);
                              zmm7 = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&im[i+16],zmm7);
                              zmm8 = _mm512_load_ps(&rho[i+32]);
                              zmm9 = _mm512_load_ps(&tht[i+32]);
                              zmm10= _mm512_mul_ps(zmm8,xcosf(zmm9));
                              _mm512_store_ps(&re[i+32],zmm10);
                              zmm11= _mm512_mul_ps(zmm8,xsinf(zmm9));
                              _mm512_store_ps(&im[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&rho[i+48]);
                              zmm13 = _mm512_load_ps(&tht[i+48]);
                              zmm14 = _mm512_mul_ps(zmm12,xcosf(zmm13));
                              _mm512_store_ps(&re[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(zmm12,xsinf(zmm13));
                              _mm512_store_ps(&im[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&rho[i+64]);
                              zmm17 = _mm512_load_ps(&tht[i+64]);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&re[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&im[i+64],zmm19);

                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0 = _mm512_load_ps(&rho[i+0]);
                              zmm1 = _mm512_load_ps(&tht[i+0]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                              _mm512_store_ps(&re[i+0],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                              _mm512_store_ps(&im[i+0],zmm3);
                              zmm4 = _mm512_load_ps(&rho[i+16]);
                              zmm5 = _mm512_load_ps(&tht[i+16]);
                              zmm6 = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&re[i+16],zmm6);
                              zmm7 = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&im[i+16],zmm7);
                              zmm8 = _mm512_load_ps(&rho[i+32]);
                              zmm9 = _mm512_load_ps(&tht[i+32]);
                              zmm10= _mm512_mul_ps(zmm8,xcosf(zmm9));
                              _mm512_store_ps(&re[i+32],zmm10);
                              zmm11= _mm512_mul_ps(zmm8,xsinf(zmm9));
                              _mm512_store_ps(&im[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&rho[i+48]);
                              zmm13 = _mm512_load_ps(&tht[i+48]);
                              zmm14 = _mm512_mul_ps(zmm12,xcosf(zmm13));
                              _mm512_store_ps(&re[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(zmm12,xsinf(zmm13));
                              _mm512_store_ps(&im[i+48],zmm15);

                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0 = _mm512_load_ps(&rho[i+0]);
                              zmm1 = _mm512_load_ps(&tht[i+0]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                              _mm512_store_ps(&re[i+0],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                              _mm512_store_ps(&im[i+0],zmm3);
                              zmm4 = _mm512_load_ps(&rho[i+16]);
                              zmm5 = _mm512_load_ps(&tht[i+16]);
                              zmm6 = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&re[i+16],zmm6);
                              zmm7 = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&im[i+16],zmm7);

                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0 = _mm512_load_ps(&rho[i+0]);
                              zmm1 = _mm512_load_ps(&tht[i+0]);
                              zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                              _mm512_store_ps(&re[i+0],zmm2);
                              zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                              _mm512_store_ps(&im[i+0],zmm3);

                         }

                          for(; (i+0) < n; i += 1) {

                               const float xr = rho[i];
                               const float xt = tht[i];
                               const float xre = xr*ceph_cosf(xt);
                               re[i] = xre;
                               const float xim = xr*ceph_sinf(xt);
                               im[i] = xim;
                         }
                         
                 }


                 

       } // math

} // gms











#endif /*__GMS_CPOLAR_ZMM16R4_HPP__*/
