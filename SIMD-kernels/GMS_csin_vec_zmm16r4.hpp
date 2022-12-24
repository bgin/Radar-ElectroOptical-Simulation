

#ifndef __GMS_CSIN_VEC_ZMM16R4_HPP__
#define __GMS_CSIN_VEC_ZMM16R4_HPP__ 231220221547

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

    const unsigned int GMS_CSIN_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CSIN_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CSIN_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CSIN_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CSIN_VEC_ZMM16R4_MAJOR+
      100U*GMS_CSIN_VEC_ZMM16R4_MINOR+
      10U*GMS_CSIN_VEC_ZMM16R4_MICRO;
    const char * const GMS_CSIN_VEC_ZMM16R4_CREATION_DATE = "23-12-2022 15:47 AM +00200 (FRI 23 DEC 2022 GMT+2)";
    const char * const GMS_CSIN_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CSIN_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CSIN_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex sine function."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"



namespace  gms {


      namespace math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csinv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict csre,
                                                   float * __restrict csim,
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
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(xsinf(zmm20),xcoshf(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(xcosf(zmm20),xsinhf(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                              zmm24 = _mm512_loadu_ps(&xre[i+96]);
                              zmm25 = _mm512_loadu_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(xsinf(zmm24),xcoshf(zmm25));
                              _mm512_storeu_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(xcosf(zmm24),xsinhf(zmm25));
                              _mm512_storeu_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_loadu_ps(&xre[i+112]);
                              zmm29 = _mm512_loadu_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(xsinf(zmm28),xcoshf(zmm29));
                              _mm512_storeu_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(xcosf(zmm28),xsinhf(zmm29));
                              _mm512_storeu_ps(&csim[i+112],zmm31);
                              _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                              zmm0 = _mm512_loadu_ps(&xre[i+128]);
                              zmm1 = _mm512_loadu_ps(&xim[i+128]);
                              zmm2 = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+128],zmm2);
                              zmm3 = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+128],zmm3);
                              zmm4 = _mm512_loadu_ps(&xre[i+144]);
                              zmm5 = _mm512_loadu_ps(&xim[i+144]);
                              zmm6 = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+144],zmm6);
                              zmm7 = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+144],zmm7);
                         }

                          for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(xsinf(zmm20),xcoshf(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(xcosf(zmm20),xsinhf(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
                              zmm24 = _mm512_loadu_ps(&xre[i+96]);
                              zmm25 = _mm512_loadu_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(xsinf(zmm24),xcoshf(zmm25));
                              _mm512_storeu_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(xcosf(zmm24),xsinhf(zmm25));
                              _mm512_storeu_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_loadu_ps(&xre[i+112]);
                              zmm29 = _mm512_loadu_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(xsinf(zmm28),xcoshf(zmm29));
                              _mm512_storeu_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(xcosf(zmm28),xsinhf(zmm29));
                              _mm512_storeu_ps(&csim[i+112],zmm31);
                         }

                          for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = ceph_sinf(z0)*ceph_coshf(z1);
                              csre[i]        = z2;
                              const float z3 = ceph_cosf(z0)*ceph_sinhf(z1);
                              csim[i]        = z3;
                         }
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csinv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) csre,
                                                   float * __restrict __ATTR_ALIGN__(64) csim,
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
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+80]);
                              zmm21 = _mm512_load_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(xsinf(zmm20),xcoshf(zmm21));
                              _mm512_store_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(xcosf(zmm20),xsinhf(zmm21));
                              _mm512_storu_ps(&csim[i+80],zmm23);
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                              zmm24 = _mm512_load_ps(&xre[i+96]);
                              zmm25 = _mm512_load_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(xsinf(zmm24),xcoshf(zmm25));
                              _mm512_store_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(xcosf(zmm24),xsinhf(zmm25));
                              _mm512_store_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_load_ps(&xre[i+112]);
                              zmm29 = _mm512_load_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(xsinf(zmm28),xcoshf(zmm29));
                              _mm512_store_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(xcosf(zmm28),xsinhf(zmm29));
                              _mm512_store_ps(&csim[i+112],zmm31);
                              _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                              zmm0 = _mm512_load_ps(&xre[i+128]);
                              zmm1 = _mm512_load_ps(&xim[i+128]);
                              zmm2 = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+128],zmm2);
                              zmm3 = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+128],zmm3);
                              zmm4 = _mm512_load_ps(&xre[i+144]);
                              zmm5 = _mm512_load_ps(&xim[i+144]);
                              zmm6 = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_store_ps(&csre[i+144],zmm6);
                              zmm7 = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_store_ps(&csim[i+144],zmm7);
                         }

                          for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+80]);
                              zmm21 = _mm512_load_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(xsinf(zmm20),xcoshf(zmm21));
                              _mm512_store_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(xcosf(zmm20),xsinhf(zmm21));
                              _mm512_store_ps(&csim[i+80],zmm23);
                              zmm24 = _mm512_load_ps(&xre[i+96]);
                              zmm25 = _mm512_load_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(xsinf(zmm24),xcoshf(zmm25));
                              _mm512_store_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(xcosf(zmm24),xsinhf(zmm25));
                              _mm512_store_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_load_ps(&xre[i+112]);
                              zmm29 = _mm512_load_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(xsinf(zmm28),xcoshf(zmm29));
                              _mm512_store_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(xcosf(zmm28),xsinhf(zmm29));
                              _mm512_store_ps(&csim[i+112],zmm31);
                         }

                          for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_load_ps(&xre[i+64]);
                              zmm17 = _mm512_load_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_store_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_store_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+32]);
                              zmm9  = _mm512_load_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_store_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_store_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_load_ps(&xre[i+48]);
                              zmm13 = _mm512_load_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_store_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_store_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_load_ps(&xre[i+16]);
                              zmm5  = _mm512_load_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_store_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_store_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_load_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_store_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_store_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = ceph_sinf(z0)*ceph_coshf(z1);
                              csre[i]        = z2;
                              const float z3 = ceph_cosf(z0)*ceph_sinhf(z1);
                              csim[i]        = z3;
                         }
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csinv_zmm16r4_unroll_8x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict csre,
                                                   float * __restrict csim,
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
 
                          for(i = 0; (i+127) < n; i += 128) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+80]);
                              zmm21 = _mm512_loadu_ps(&xim[i+80]);
                              zmm22 = _mm512_mul_ps(xsinf(zmm20),xcoshf(zmm21));
                              _mm512_storeu_ps(&csre[i+80],zmm22);
                              zmm23 = _mm512_mul_ps(xcosf(zmm20),xsinhf(zmm21));
                              _mm512_storeu_ps(&csim[i+80],zmm23);
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                              zmm24 = _mm512_loadu_ps(&xre[i+96]);
                              zmm25 = _mm512_loadu_ps(&xim[i+96]);
                              zmm26 = _mm512_mul_ps(xsinf(zmm24),xcoshf(zmm25));
                              _mm512_storeu_ps(&csre[i+96],zmm26);
                              zmm27 = _mm512_mul_ps(xcosf(zmm24),xsinhf(zmm25));
                              _mm512_storeu_ps(&csim[i+96],zmm27);
                              zmm28 = _mm512_loadu_ps(&xre[i+112]);
                              zmm29 = _mm512_loadu_ps(&xim[i+112]);
                              zmm30 = _mm512_mul_ps(xsinf(zmm28),xcoshf(zmm29));
                              _mm512_storeu_ps(&csre[i+112],zmm30);
                              zmm31 = _mm512_mul_ps(xcosf(zmm28),xsinhf(zmm29));
                              _mm512_storeu_ps(&csim[i+112],zmm31);
                             
                         }

                         
                          for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                              zmm16 = _mm512_loadu_ps(&xre[i+64]);
                              zmm17 = _mm512_loadu_ps(&xim[i+64]);
                              zmm18 = _mm512_mul_ps(xsinf(zmm16),xcoshf(zmm17));
                              _mm512_storeu_ps(&csre[i+64],zmm18);
                              zmm19 = _mm512_mul_ps(xcosf(zmm16),xsinhf(zmm17));
                              _mm512_storeu_ps(&csim[i+64],zmm19);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+32]);
                              zmm9  = _mm512_loadu_ps(&xim[i+32]);
                              zmm10 = _mm512_mul_ps(xsinf(zmm8),xcoshf(zmm9));
                              _mm512_storeu_ps(&csre[i+32],zmm10);
                              zmm11 = _mm512_mul_ps(xcosf(zmm8),xsinhf(zmm9));
                              _mm512_storeu_ps(&csim[i+32],zmm11);
                              zmm12 = _mm512_loadu_ps(&xre[i+48]);
                              zmm13 = _mm512_loadu_ps(&xim[i+48]);
                              zmm14 = _mm512_mul_ps(xsinf(zmm12),xcoshf(zmm13));
                              _mm512_storeu_ps(&csre[i+48],zmm14);
                              zmm15 = _mm512_mul_ps(xcosf(zmm12),xsinhf(zmm13));
                              _mm512_storeu_ps(&csim[i+48],zmm15);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                              zmm4  = _mm512_loadu_ps(&xre[i+16]);
                              zmm5  = _mm512_loadu_ps(&xim[i+16]);
                              zmm6  = _mm512_mul_ps(xsinf(zmm4),xcoshf(zmm5));
                              _mm512_storeu_ps(&csre[i+16],zmm6);
                              zmm7  = _mm512_mul_ps(xcosf(zmm4),xsinhf(zmm5));
                              _mm512_storeu_ps(&csim[i+16],zmm7);
                         }

                          for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_loadu_ps(&xim[i+0]);
                              zmm2  = _mm512_mul_ps(xsinf(zmm0),xcoshf(zmm1));
                              _mm512_storeu_ps(&csre[i+0],zmm2);
                              zmm3  = _mm512_mul_ps(xcosf(zmm0),xsinhf(zmm1));
                              _mm512_storeu_ps(&csim[i+0],zmm3);
                         }

                          for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = xim[i];
                              const float z2 = ceph_sinf(z0)*ceph_coshf(z1);
                              csre[i]        = z2;
                              const float z3 = ceph_cosf(z0)*ceph_sinhf(z1);
                              csim[i]        = z3;
                         }
                }




     }

} // gms





#endif /*__GMS_CSIN_VEC_ZMM16R4_HPP__*/
