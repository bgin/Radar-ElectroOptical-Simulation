




#include <immintrin.h>
#include "GMS_cconj_zmm16r4.h"




                   void gms::math::cconjv_zmm16r4_unroll_16x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                     
                      
                      int32_t i;
                      if(__builtin_expect(0==n,0))
                      {
                         return;
                      }
                      else if(n>0 && n<16)
                      {
                          for(i = 0;i < n; i++) {
                           const float z0 = xim[i];
                           cxim[i]        = -1.0f*z0;
                          }
                          return;
                      }
                      else if(n<=32) 
                      {
                          const __m256 CN1   = _mm256_set1_ps(-0x1p+0);
                          _mm_prefetch((consrt char*)&xim[0],_MM_HINT_T0);
                          __m256 ymm0 = _mm256_loadu_ps(&xim[0]);
                          _mm256_storeu_ps(&cxim[0], _mm256_mul_ps(CN1,ymm0));
                          __m256 ymm1 = _mm256_loadu_ps(&xim[8]);
                          _mm256_storeu_ps(&cxim[8], _mm256_mul_ps(CN1,ymm1));
                          __m256 ymm2 = _mm256_loadu_ps(&xim[16]);
                          _mm256_storeu_ps(&cxim[16], _mm256_mul_ps(CN1,ymm2));
                          __m256 ymm3 = _mm256_loadu_ps(&xim[24]);
                          _mm256_storeu_ps(&cxim[24], _mm256_mul_ps(CN1,ymm3));
                          return;
                      }
                      else if(n<=64)
                      {
                           const __m256 CN1   = _mm256_set1_ps(-0x1p+0);
                            _mm_prefetch((consrt char*)&xim[0],_MM_HINT_T0);
                          __m256 ymm0 = _mm256_loadu_ps(&xim[0]);
                          _mm256_storeu_ps(&cxim[0], _mm256_mul_ps(CN1,ymm0));
                          __m256 ymm1 = _mm256_loadu_ps(&xim[8]);
                          _mm256_storeu_ps(&cxim[8], _mm256_mul_ps(CN1,ymm1));
                          __m256 ymm2 = _mm256_loadu_ps(&xim[16]);
                          _mm256_storeu_ps(&cxim[16], _mm256_mul_ps(CN1,ymm2));
                          __m256 ymm3 = _mm256_loadu_ps(&xim[24]);
                           _mm_prefetch((consrt char*)&xim[32],_MM_HINT_T0);
                          _mm256_storeu_ps(&cxim[24], _mm256_mul_ps(CN1,ymm3));
                          __m256 ymm4 = _mm256_loadu_ps(&xim[32]);
                          _mm256_storeu_ps(&cxim[32], _mm256_mul_ps(CN1,ymm4));
                          __m256 ymm5 = _mm256_loadu_ps(&xim[40]);
                          _mm256_storeu_ps(&cxim[40], _mm256_mul_ps(CN1,ymm5));
                          __m256 ymm6 = _mm256_loadu_ps(&xim[48]);
                          _mm256_storeu_ps(&cxim[48], _mm256_mul_ps(CN1,ymm6));
                          __m256 ymm7 = _mm256_loadu_ps(&xim[56]);
                          _mm256_storeu_ps(&cxim[56], _mm256_mul_ps(CN1,ymm7));
                          return;
                      }
                      else if(n>64) {
                           __m512 zmm0,zmm1,zmm2,zmm3;
                           __m512 zmm4,zmm5,zmm6,zmm7;
                           __m512 zmm8,zmm9,zmm10,zmm11;
                           __m512 zmm12,zmm13,zmm14,zmm15;
                           const  __m512 none = _mm512_set1_ps(-0x1p+0);
                          for(i = 0; (i+255) < n; i += 256) {
                          _mm_prefetch((const char*)&xim[i+256], _MM_HINT_T0);
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
                          zmm12 = _mm512_loadu_ps(&xim[i+192]);
                          _mm512_storeu_ps(&cxim[i+192], _mm512_mul_ps(none,zmm12));
                          zmm13 = _mm512_loadu_ps(&xim[i+208]);
                          _mm512_storeu_ps(&cxim[i+208], _mm512_mul_ps(none,zmm13));
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
               }


                 
                   void gms::math::cconjv_zmm16r4_unroll_16x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                       __m512 zmm0,zmm1,zmm2,zmm3;
                       __m512 zmm4,zmm5,zmm6,zmm7;
                       __m512 zmm8,zmm9,zmm10,zmm11;
                       __m512 zmm12,zmm13,zmm14,zmm15;
                      const  __m512 none = _mm512_set1_ps(-1.0f);
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


                  
                   void gms::math::cconjv_zmm16r4_unroll_10x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                       __m512 zmm0,zmm1,zmm2,zmm3;
                       __m512 zmm4,zmm5,zmm6,zmm7;
                       __m512 zmm8,zmm9;
                      const  __m512 none = _mm512_set1_ps(-1.0f);
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


                  
                   void gms::math::cconjv_zmm16r4_unroll_10x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                       __m512 zmm0,zmm1,zmm2,zmm3;
                       __m512 zmm4,zmm5,zmm6,zmm7;
                       __m512 zmm8,zmm9;
                      const  __m512 none = _mm512_set1_ps(-1.0f);
                      int32_t i;

                      for(i = 0; (i+159) < n; i += 160) {
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


                  
                   void gms::math::cconjv_zmm16r4_unroll_6x_u(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                       __m512 zmm0,zmm1,zmm2,zmm3;
                       __m512 zmm4,zmm5;
                      const  __m512 none = _mm512_set1_ps(-1.0f);
                      int32_t i;

                      for(i = 0; (i+95) < n; i += 96) {
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


                 
                   void gms::math::cconjv_zmm16r4_unroll_6x_a(const float * __restrict xim,
                                                    float * __restrict cxim,
                                                    const int32_t n) {

                      if(__builtin_expect(0==n,0)) {return;}
                       __m512 zmm0,zmm1,zmm2,zmm3;
                       __m512 zmm4,zmm5;
                      const  __m512 none = _mm512_set1_ps(-1.0f);
                      int32_t i;

                      for(i = 0; (i+95) < n; i += 96) {
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




                  
                   

        
