
#ifndef __GMS_CCONJ_VEC_ZMM16R4_HPP__
#define __GMS_CCONJ_VEC_ZMM16R4_HPP__ 281220221036


namespace file_version {

    const unsigned int GMS_CCONJ_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CCONJ_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CCONJ_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CCONJ_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CCONJ_VEC_ZMM16R4_MAJOR+
      100U*GMS_CCONJ_VEC_ZMM16R4_MINOR+
      10U*GMS_CCONJ_VEC_ZMM16R4_MICRO;
    const char * const GMS_CCONJ_VEC_ZMM16R4_CREATION_DATE = "28-12-2022 10:36  +00200 (WED 28 DEC 2022 GMT+2)";
    const char * const GMS_CCONJ_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CCONJ_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CCONJ_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector conjugate operations."

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
                   void cconjv_zmm16r4_unroll_16x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      register __m512 zmm4,zmm5,zmm6,zmm7;
                      register __m512 zmm8,zmm9,zmm10,zmm11;
                      register __m512 zmm12,zmm13,zmm14,zmm15;
                      const register __m512 none = _mm512_set1_ps(-1.0f);
                      int32_t i;

                      for(i = 0; (i+255) < n; i += 256) {
                          _mm_prefetch((const char*)&xim[i+32], _MM_HINT_T0);
                          zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_loadu_ps(&xim[i+16]);
                          _mm512_storeu_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          _mm_prefetch((const char*)&xim[i+64], _MM_HINT_T0);
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_storeu_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_loadu_ps(&xim[i+48]);
                          _mm512_storeu_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          _mm_prefetch((const char*)&xim[i+96], _MM_HINT_T0);
                          zmm4 = _mm512_loadu_ps(&xim[i+64]);
                          _mm512_storeu_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_loadu_ps(&xim[i+80]);
                          _mm512_storeu_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          _mm_prefetch((const char*)&xim[i+128], _MM_HINT_T0);
                          zmm6 = _mm512_loadu_ps(&xim[i+96]);
                          _mm512_storeu_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_loadu_ps(&xim[i+112]);
                          _mm512_storeu_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                          _mm_prefetch((const char*)&xim[i+160], _MM_HINT_T0);
                          zmm8 = _mm512_loadu_ps(&xim[i+128]);
                          _mm512_storeu_ps(&cxim[i+128], _mm512_mul_ps(none,zmm8));
                          zmm9 = _mm512_loadu_ps(&xim[i+144]);
                          _mm512_storeu_ps(&cxim[i+144], _mm512_mul_ps(none,zmm9));
                          _mm_prefetch((const char*)&xim[i+192], _MM_HINT_T0);
                          zmm10 = _mm512_loadu_ps(&xim[i+160]);
                          _mm512_storeu_ps(&cxim[i+160], _mm512_mul_ps(none,zmm10));
                          zmm11 = _mm512_loadu_ps(&xim[i+176]);
                          _mm512_storeu_ps(&cxim[i+176], _mm512_mul_ps(none,zmm11));
                          _mm_prefetch((const char*)&xim[i+224], _MM_HINT_T0);
                          zmm12 = _mm512_loadu_ps(&xim[i+192]);
                          _mm512_storeu_ps(&cxim[i+192], _mm512_mul_ps(none,zmm12));
                          zmm13 = _mm512_loadu_ps(&xim[i+208]);
                          _mm512_storeu_ps(&cxim[i+208], _mm512_mul_ps(none,zmm13));
                          _mm_prefetch((const char*)&xim[i+256], _MM_HINT_T0);
                          zmm14 = _mm512_loadu_ps(&xim[i+224]);
                          _mm512_storeu_ps(&cxim[i+224], _mm512_mul_ps(none,zmm14));
                          zmm15 = _mm512_loadu_ps(&xim[i+240]);
                          _mm512_storeu_ps(&cxim[i+240], _mm512_mul_ps(none,zmm15));
                     }

                       for(; (i+191) < n; i += 192) {
                            zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_loadu_ps(&xim[i+16]);
                          _mm512_storeu_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_storeu_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_loadu_ps(&xim[i+48]);
                          _mm512_storeu_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_loadu_ps(&xim[i+64]);
                          _mm512_storeu_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_loadu_ps(&xim[i+80]);
                          _mm512_storeu_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          zmm6 = _mm512_loadu_ps(&xim[i+96]);
                          _mm512_storeu_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_loadu_ps(&xim[i+112]);
                          _mm512_storeu_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                          zmm8 = _mm512_loadu_ps(&xim[i+128]);
                          _mm512_storeu_ps(&cxim[i+128], _mm512_mul_ps(none,zmm8));
                          zmm9 = _mm512_loadu_ps(&xim[i+144]);
                          _mm512_storeu_ps(&cxim[i+144], _mm512_mul_ps(none,zmm9));
                          zmm10 = _mm512_loadu_ps(&xim[i+160]);
                          _mm512_storeu_ps(&cxim[i+160], _mm512_mul_ps(none,zmm10));
                          zmm11 = _mm512_loadu_ps(&xim[i+176]);
                          _mm512_storeu_ps(&cxim[i+176], _mm512_mul_ps(none,zmm11));
                     }

                       for(; (i+127) < n; i += 128) {
                           zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_loadu_ps(&xim[i+16]);
                          _mm512_storeu_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_storeu_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_loadu_ps(&xim[i+48]);
                          _mm512_storeu_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_loadu_ps(&xim[i+64]);
                          _mm512_storeu_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_loadu_ps(&xim[i+80]);
                          _mm512_storeu_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          zmm6 = _mm512_loadu_ps(&xim[i+96]);
                          _mm512_storeu_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_loadu_ps(&xim[i+112]);
                          _mm512_storeu_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                     }

                       for(; (i+79) < n; i += 80) {
                           zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_loadu_ps(&xim[i+16]);
                          _mm512_storeu_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_storeu_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_loadu_ps(&xim[i+48]);
                          _mm512_storeu_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_loadu_ps(&xim[i+64]);
                          _mm512_storeu_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                     }

                       for(; (i+63) < n; i += 64) {
                          zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_loadu_ps(&xim[i+16]);
                          _mm512_storeu_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_storeu_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_loadu_ps(&xim[i+48]);
                          _mm512_storeu_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3)); 
                     }

                       for(; (i+31) < n; i += 32) {
                          zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_loadu_ps(&xim[i+16]);
                          _mm512_storeu_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                     }

                       for(; (i+15) < n; i += 16) {
                           zmm0 = _mm512_loadu_ps(&xim[i+0]);
                          _mm512_storeu_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                     }

                       for(; (i+0) < n; i += 1) {
                           const float z0 = xim[i];
                           cxim[i]        = -1.0f*z0;
                     }
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cconjv_zmm16r4_unroll_16x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      register __m512 zmm4,zmm5,zmm6,zmm7;
                      register __m512 zmm8,zmm9,zmm10,zmm11;
                      register __m512 zmm12,zmm13,zmm14,zmm15;
                      const register __m512 none = _mm512_set1_ps(-1.0f);
                      int32_t i;

                      for(i = 0; (i+255) < n; i += 256) {
                          _mm_prefetch((const char*)&xim[i+32], _MM_HINT_T0);
                          zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          _mm_prefetch((const char*)&xim[i+64], _MM_HINT_T0);
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          _mm_prefetch((const char*)&xim[i+96], _MM_HINT_T0);
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_load_ps(&xim[i+80]);
                          _mm512_store_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          _mm_prefetch((const char*)&xim[i+128], _MM_HINT_T0);
                          zmm6 = _mm512_load_ps(&xim[i+96]);
                          _mm512_storeu_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_load_ps(&xim[i+112]);
                          _mm512_store_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                          _mm_prefetch((const char*)&xim[i+160], _MM_HINT_T0);
                          zmm8 = _mm512_load_ps(&xim[i+128]);
                          _mm512_store_ps(&cxim[i+128], _mm512_mul_ps(none,zmm8));
                          zmm9 = _mm512_load_ps(&xim[i+144]);
                          _mm512_store_ps(&cxim[i+144], _mm512_mul_ps(none,zmm9));
                          _mm_prefetch((const char*)&xim[i+192], _MM_HINT_T0);
                          zmm10 = _mm512_load_ps(&xim[i+160]);
                          _mm512_store_ps(&cxim[i+160], _mm512_mul_ps(none,zmm10));
                          zmm11 = _mm512_load_ps(&xim[i+176]);
                          _mm512_store_ps(&cxim[i+176], _mm512_mul_ps(none,zmm11));
                          _mm_prefetch((const char*)&xim[i+224], _MM_HINT_T0);
                          zmm12 = _mm512_load_ps(&xim[i+192]);
                          _mm512_store_ps(&cxim[i+192], _mm512_mul_ps(none,zmm12));
                          zmm13 = _mm512_load_ps(&xim[i+208]);
                          _mm512_store_ps(&cxim[i+208], _mm512_mul_ps(none,zmm13));
                          _mm_prefetch((const char*)&xim[i+256], _MM_HINT_T0);
                          zmm14 = _mm512_load_ps(&xim[i+224]);
                          _mm512_store_ps(&cxim[i+224], _mm512_mul_ps(none,zmm14));
                          zmm15 = _mm512_load_ps(&xim[i+240]);
                          _mm512_store_ps(&cxim[i+240], _mm512_mul_ps(none,zmm15));
                     }

                       for(; (i+191) < n; i += 192) {
                            zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_load_ps(&xim[i+80]);
                          _mm512_store_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          zmm6 = _mm512_load_ps(&xim[i+96]);
                          _mm512_store_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_load_ps(&xim[i+112]);
                          _mm512_store_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                          zmm8 = _mm512_load_ps(&xim[i+128]);
                          _mm512_store_ps(&cxim[i+128], _mm512_mul_ps(none,zmm8));
                          zmm9 = _mm512_load_ps(&xim[i+144]);
                          _mm512_store_ps(&cxim[i+144], _mm512_mul_ps(none,zmm9));
                          zmm10 = _mm512_load_ps(&xim[i+160]);
                          _mm512_store_ps(&cxim[i+160], _mm512_mul_ps(none,zmm10));
                          zmm11 = _mm512_load_ps(&xim[i+176]);
                          _mm512_store_ps(&cxim[i+176], _mm512_mul_ps(none,zmm11));
                     }

                       for(; (i+127) < n; i += 128) {
                           zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_load_ps(&xim[i+80]);
                          _mm512_store_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          zmm6 = _mm512_load_ps(&xim[i+96]);
                          _mm512_store_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_load_ps(&xim[i+112]);
                          _mm512_store_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                     }

                       for(; (i+79) < n; i += 80) {
                           zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                     }

                       for(; (i+63) < n; i += 64) {
                          zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3)); 
                     }

                       for(; (i+31) < n; i += 32) {
                          zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                     }

                       for(; (i+15) < n; i += 16) {
                           zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                     }

                       for(; (i+0) < n; i += 1) {
                           const float z0 = xim[i];
                           cxim[i]        = -1.0f*z0;
                     }
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void cconjv_zmm16r4_unroll_10x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      register __m512 zmm4,zmm5,zmm6,zmm7;
                      register __m512 zmm8,zmm9;
                      const register __m512 none = _mm512_set1_ps(-1.0f);
                      int32_t i;

                      for(i = 0; (i+159) < n; i += 160) {
                          _mm_prefetch((const char*)&xim[i+32], _MM_HINT_T0);
                          zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          _mm_prefetch((const char*)&xim[i+64], _MM_HINT_T0);
                          zmm2 = _mm512_loadu_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          _mm_prefetch((const char*)&xim[i+96], _MM_HINT_T0);
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_load_ps(&xim[i+80]);
                          _mm512_store_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          _mm_prefetch((const char*)&xim[i+128], _MM_HINT_T0);
                          zmm6 = _mm512_load_ps(&xim[i+96]);
                          _mm512_storeu_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_load_ps(&xim[i+112]);
                          _mm512_store_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                          _mm_prefetch((const char*)&xim[i+160], _MM_HINT_T0);
                          zmm8 = _mm512_load_ps(&xim[i+128]);
                          _mm512_store_ps(&cxim[i+128], _mm512_mul_ps(none,zmm8));
                          zmm9 = _mm512_load_ps(&xim[i+144]);
                          _mm512_store_ps(&cxim[i+144], _mm512_mul_ps(none,zmm9));
                          
                     }

                      
                       for(; (i+127) < n; i += 128) {
                           zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                          zmm5 = _mm512_load_ps(&xim[i+80]);
                          _mm512_store_ps(&cxim[i+80], _mm512_mul_ps(none,zmm5));
                          zmm6 = _mm512_load_ps(&xim[i+96]);
                          _mm512_store_ps(&cxim[i+96], _mm512_mul_ps(none,zmm6));
                          zmm7 = _mm512_load_ps(&xim[i+112]);
                          _mm512_store_ps(&cxim[i+112], _mm512_mul_ps(none,zmm7));
                     }

                       for(; (i+79) < n; i += 80) {
                           zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3));
                          zmm4 = _mm512_load_ps(&xim[i+64]);
                          _mm512_store_ps(&cxim[i+64], _mm512_mul_ps(none,zmm4));
                     }

                       for(; (i+63) < n; i += 64) {
                          zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                          zmm2 = _mm512_load_ps(&xim[i+32]);
                          _mm512_store_ps(&cxim[i+32], _mm512_mul_ps(none,zmm2));
                          zmm3 = _mm512_load_ps(&xim[i+48]);
                          _mm512_store_ps(&cxim[i+48], _mm512_mul_ps(none,zmm3)); 
                     }

                       for(; (i+31) < n; i += 32) {
                          zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                          zmm1 = _mm512_load_ps(&xim[i+16]);
                          _mm512_store_ps(&cxim[i+16], _mm512_mul_ps(none,zmm1));
                     }

                       for(; (i+15) < n; i += 16) {
                           zmm0 = _mm512_load_ps(&xim[i+0]);
                          _mm512_store_ps(&cxim[i+0], _mm512_mul_ps(none,zmm0));
                     }

                       for(; (i+0) < n; i += 1) {
                           const float z0 = xim[i];
                           cxim[i]        = -1.0f*z0;
                     }
               }



        } // math


} // gms















#endif /*__GMS_CCONJ_VEC_ZMM16R4_HPP__*/
