

#ifndef __GMS_CSQRT_VEC_ZMM16R4_HPP__
#define __GMS_CSQRT_VEC_ZMM16R4_HPP__ 271220221650


namespace file_version {

    const unsigned int GMS_CSQRT_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CSQRT_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CSQRT_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CSQRT_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CSQRT_VEC_ZMM16R4_MAJOR+
      100U*GMS_CSQRT_VEC_ZMM16R4_MINOR+
      10U*GMS_CSQRT_VEC_ZMM16R4_MICRO;
    const char * const GMS_CSQRT_VEC_ZMM16R4_CREATION_DATE = "27-12-2022 16:50  +00200 (TUE 27 DEC 2022 GMT+2)";
    const char * const GMS_CSQRT_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CSQRT_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CSQRT_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector sqrt operations."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_cephes.h"
#include "GMS_cabs_vec_zmm16r4.hpp"

namespace gms {

        namespace math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void csqrtv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                    const float * __restrict yim,
                                                    float * __restrict wrkc,
                                                    float * __restrict csqr,
                                                    float * __restrict csqi,
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
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cabsv_zmm16r4_unroll_10x_u(xre,xim,wrkc,n);
                        for(i = 0; (i+159) < n; i += 160) {
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrkc[i+32],_MM_HINT_T0);
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&wrkc[i+0]);
                            zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+0],_mm512_sqrt_ps(zmm2));
                            zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+0],_mm512_sqrt_ps(zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&wrkc[i+16]);
                            zmm6  = _mm512_mul_ps(half,_mm512_add_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqr[i+16],_mm512_sqrt_ps(zmm6));
                            zmm7  = _mm512_mul_ps(half,_mm512_sub_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqi[i+16],_mm512_sqrt_ps(zmm7));
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrkc[i+64],_MM_HINT_T0);
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&wrkc[i+32]);
                            zmm10 = _mm512_mul_ps(half,_mm512_add_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqr[i+32],_mm512_sqrt_ps(zmm10));
                            zmm11 = _mm512_mul_ps(half,_mm512_sub_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqi[i+32],_mm512_sqrt_ps(zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&wrkc[i+48]);
                            zmm14 = _mm512_mul_ps(half,_mm512_add_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqr[i+48],_mm512_sqrt_ps(zmm14));
                            zmm15 = _mm512_mul_ps(half,_mm512_sub_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqi[i+48],_mm512_sqrt_ps(zmm15));
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrkc[i+96],_MM_HINT_T0);
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&wrkc[i+64]);
                            zmm18 = _mm512_mul_ps(half,_mm512_add_ps(zmm17,zmm16));
                            _mm512_storeu_ps(&csqr[i+64],_mm512_sqrt_ps(zmm18));
                            zmm19 = _mm512_mul_ps(half,_mm512_sub_ps(zmm17,zmm16));
                            _mm512_storeu_ps(&csqi[i+64],_mm512_sqrt_ps(zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&wrkc[i+80]);
                            zmm22 = _mm512_mul_ps(half,_mm512_add_ps(zmm21,zmm20));
                            _mm512_storeu_ps(&csqr[i+80],_mm512_sqrt_ps(zmm22));
                            zmm23 = _mm512_mul_ps(half,_mm512_sub_ps(zmm21,zmm20));
                            _mm512_storeu_ps(&csqi[i+80],_mm512_sqrt_ps(zmm23));
                            _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrkc[i+128],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&xre[i+96]);
                            zmm25 = _mm512_loadu_ps(&wrkc[i+96]);
                            zmm26 = _mm512_mul_ps(half,_mm512_add_ps(zmm25,zmm24));
                            _mm512_storeu_ps(&csqr[i+96],_mm512_sqrt_ps(zmm26));
                            zmm27 = _mm512_mul_ps(half,_mm512_sub_ps(zmm25,zmm24));
                            _mm512_storeu_ps(&csqi[i+96],_mm512_sqrt_ps(zmm27));
                            zmm28 = _mm512_loadu_ps(&xre[i+112]);
                            zmm29 = _mm512_loadu_ps(&wrkc[i+112]);
                            zmm30 = _mm512_mul_ps(half,_mm512_add_ps(zmm29,zmm28));
                            _mm512_storeu_ps(&csqr[i+112],_mm512_sqrt_ps(zmm30));
                            zmm31 = _mm512_mul_ps(half,_mm512_sub_ps(zmm29,zmm28));
                            _mm512_storeu_ps(&csqi[i+112],_mm512_sqrt_ps(zmm31));
                            _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrkc[i+160],_MM_HINT_T0);
                            zmm0 = _mm512_loadu_ps(&xre[i+128]);
                            zmm1 = _mm512_loadu_ps(&wrkc[i+128]);
                            zmm2 = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+128],_mm512_sqrt_ps(zmm2));
                            zmm3 = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+128],_mm512_sqrt_ps(zmm3));
                            zmm4 = _mm512_loadu_ps(&xre[i+144]);
                            zmm5 = _mm512_loadu_ps(&wrkc[i+144]);
                            zmm6 = _mm512_mul_ps(half,_mm512_add_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqr[i+144],_mm512_sqrt_ps(zmm6));
                            zmm7 = _mm512_mul_ps(half,_mm512_sub_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqi[i+144],_mm512_sqrt_ps(zmm7));
                      }

                       for(; (i+95) < n; i += 96) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&wrkc[i+0]);
                            zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+0],_mm512_sqrt_ps(zmm2));
                            zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+0],_mm512_sqrt_ps(zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&wrkc[i+16]);
                            zmm6  = _mm512_mul_ps(half,_mm512_add_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqr[i+16],_mm512_sqrt_ps(zmm6));
                            zmm7  = _mm512_mul_ps(half,_mm512_sub_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqi[i+16],_mm512_sqrt_ps(zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&wrkc[i+32]);
                            zmm10 = _mm512_mul_ps(half,_mm512_add_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqr[i+32],_mm512_sqrt_ps(zmm10));
                            zmm11 = _mm512_mul_ps(half,_mm512_sub_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqi[i+32],_mm512_sqrt_ps(zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&wrkc[i+48]);
                            zmm14 = _mm512_mul_ps(half,_mm512_add_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqr[i+48],_mm512_sqrt_ps(zmm14));
                            zmm15 = _mm512_mul_ps(half,_mm512_sub_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqi[i+48],_mm512_sqrt_ps(zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&wrkc[i+64]);
                            zmm18 = _mm512_mul_ps(half,_mm512_add_ps(zmm17,zmm16));
                            _mm512_storeu_ps(&csqr[i+64],_mm512_sqrt_ps(zmm18));
                            zmm19 = _mm512_mul_ps(half,_mm512_sub_ps(zmm17,zmm16));
                            _mm512_storeu_ps(&csqi[i+64],_mm512_sqrt_ps(zmm19));
                            zmm20 = _mm512_loadu_ps(&xre[i+80]);
                            zmm21 = _mm512_loadu_ps(&wrkc[i+80]);
                            zmm22 = _mm512_mul_ps(half,_mm512_add_ps(zmm21,zmm20));
                            _mm512_storeu_ps(&csqr[i+80],_mm512_sqrt_ps(zmm22));
                            zmm23 = _mm512_mul_ps(half,_mm512_sub_ps(zmm21,zmm20));
                            _mm512_storeu_ps(&csqi[i+80],_mm512_sqrt_ps(zmm23));

                      }

                        for(; (i+79) < n; i += 80) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&wrkc[i+0]);
                            zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+0],_mm512_sqrt_ps(zmm2));
                            zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+0],_mm512_sqrt_ps(zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&wrkc[i+16]);
                            zmm6  = _mm512_mul_ps(half,_mm512_add_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqr[i+16],_mm512_sqrt_ps(zmm6));
                            zmm7  = _mm512_mul_ps(half,_mm512_sub_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqi[i+16],_mm512_sqrt_ps(zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&wrkc[i+32]);
                            zmm10 = _mm512_mul_ps(half,_mm512_add_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqr[i+32],_mm512_sqrt_ps(zmm10));
                            zmm11 = _mm512_mul_ps(half,_mm512_sub_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqi[i+32],_mm512_sqrt_ps(zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&wrkc[i+48]);
                            zmm14 = _mm512_mul_ps(half,_mm512_add_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqr[i+48],_mm512_sqrt_ps(zmm14));
                            zmm15 = _mm512_mul_ps(half,_mm512_sub_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqi[i+48],_mm512_sqrt_ps(zmm15));
                            zmm16 = _mm512_loadu_ps(&xre[i+64]);
                            zmm17 = _mm512_loadu_ps(&wrkc[i+64]);
                            zmm18 = _mm512_mul_ps(half,_mm512_add_ps(zmm17,zmm16));
                            _mm512_storeu_ps(&csqr[i+64],_mm512_sqrt_ps(zmm18));
                            zmm19 = _mm512_mul_ps(half,_mm512_sub_ps(zmm17,zmm16));
                            _mm512_storeu_ps(&csqi[i+64],_mm512_sqrt_ps(zmm19));
                      }

                        for(; (i+63) < n; i += 64) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&wrkc[i+0]);
                            zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+0],_mm512_sqrt_ps(zmm2));
                            zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+0],_mm512_sqrt_ps(zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&wrkc[i+16]);
                            zmm6  = _mm512_mul_ps(half,_mm512_add_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqr[i+16],_mm512_sqrt_ps(zmm6));
                            zmm7  = _mm512_mul_ps(half,_mm512_sub_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqi[i+16],_mm512_sqrt_ps(zmm7));
                            zmm8  = _mm512_loadu_ps(&xre[i+32]);
                            zmm9  = _mm512_loadu_ps(&wrkc[i+32]);
                            zmm10 = _mm512_mul_ps(half,_mm512_add_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqr[i+32],_mm512_sqrt_ps(zmm10));
                            zmm11 = _mm512_mul_ps(half,_mm512_sub_ps(zmm9,zmm8));
                            _mm512_storeu_ps(&csqi[i+32],_mm512_sqrt_ps(zmm11));
                            zmm12 = _mm512_loadu_ps(&xre[i+48]);
                            zmm13 = _mm512_loadu_ps(&wrkc[i+48]);
                            zmm14 = _mm512_mul_ps(half,_mm512_add_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqr[i+48],_mm512_sqrt_ps(zmm14));
                            zmm15 = _mm512_mul_ps(half,_mm512_sub_ps(zmm13,zmm12));
                            _mm512_storeu_ps(&csqi[i+48],_mm512_sqrt_ps(zmm15));
                      }

                        for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&wrkc[i+0]);
                            zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+0],_mm512_sqrt_ps(zmm2));
                            zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+0],_mm512_sqrt_ps(zmm3));
                            zmm4  = _mm512_loadu_ps(&xre[i+16]);
                            zmm5  = _mm512_loadu_ps(&wrkc[i+16]);
                            zmm6  = _mm512_mul_ps(half,_mm512_add_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqr[i+16],_mm512_sqrt_ps(zmm6));
                            zmm7  = _mm512_mul_ps(half,_mm512_sub_ps(zmm5,zmm4));
                            _mm512_storeu_ps(&csqi[i+16],_mm512_sqrt_ps(zmm7));
                     }

                       for(; (i+15) < n; i += 16) {
                            zmm0  = _mm512_loadu_ps(&xre[i+0]);
                            zmm1  = _mm512_loadu_ps(&wrkc[i+0]);
                            zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqr[i+0],_mm512_sqrt_ps(zmm2));
                            zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                            _mm512_storeu_ps(&csqi[i+0],_mm512_sqrt_ps(zmm3));
                     }

                       for(; (i+0) < n; i += 1) {
                            const float z0 = xre[i];
                            const float z1 = wrkc[i];
                            const float z2 = 0.5f*(z1+z0);
                            csqr[i]        = ceph_sqrtf(z2);
                            const float z3 = 0.5f*(z1-z0);
                            csqi[i]        = ceph_sqrtf(z3);
                     }
                }




      } // math


} // gms
















#endif /*__GMS_CSQRT_VEC_ZMM16R4_HPP__*/
