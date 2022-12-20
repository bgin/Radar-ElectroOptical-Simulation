

#ifndef __GMS_CLOG_VEC_ZMM16R4_HPP__
#define __GMS_CLOG_VEC_ZMM16R4_HPP__ 201220220929


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

    const unsigned int GMS_CLOG_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CLOG_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CLOG_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CLOG_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CLOG_VEC_ZMM16R4_MAJOR+
      100U*GMS_CLOG_VEC_ZMM16R4_MINOR+
      10U*GMS_CLOG_VEC_ZMM16R4_MICRO;
    const char * const GMS_CLOG_VEC_ZMM16R4_CREATION_DATE = "20-12-2022 09:29 AM +00200 (TUE 20 DEC 2022 GMT+2)";
    const char * const GMS_CLOG_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CLOG_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CLOG_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex logarithm function."

}


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"
#include "GMS_carg_vec_zmm16r4.hpp"
#include "GMS_cabs_vec_zmm16r4.hpp"


namespace  gms {


         namespace math {

            
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void clogv_zmm16r4_unroll_16x_u(const float * __restrict re,
                                                   const float * __restrict im,
                                                   float * __restrict wrk1,
                                                   float * __restrict wrk2,
                                                   float * __restrict clogr,
                                                   float * __restrict clogi,
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
                          cargv_zmm16r4_unroll_16x_u(re,im,wrk1,n);
                          cabsv_zmm16r4_unroll_16x_u(re,im,wrk2,n);
                          // wrk1 -- cabs
                          // wrk2 -- carg
                        for(; (i+255) < n; i += 256) {
                            _mm_prefetch((const char*)&wrk1[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            _mm_prefetch((const char*)&wrk1[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+64],_MM_HINT_T0);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            _mm_prefetch((const char*)&wrk1[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                            _mm512_storeu_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_loadu_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_storeu_ps(&clogr[i+80],zmm16);
                            _mm512_storeu_ps(&clogi[i+80],zmm17);
                            _mm_prefetch((const char*)&wrk1[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+128],_MM_HINT_T0);
                            zmm18 = _mm512_loadu_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_loadu_ps(&wrk2[i+96]);
                            _mm512_storeu_ps(&clogr[i+96],zmm19);
                            _mm512_storeu_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_loadu_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_loadu_ps(&wrk2[i+112]);
                            _mm512_storeu_ps(&clogr[i+112],zmm22);
                            _mm512_storeu_ps(&clogi[i+112],zmm23);
                            _mm_prefetch((const char*)&wrk1[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+160],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&wrk1[i+128]);
                            zmm25 = xlogf(zmm24);
                            zmm26 = _mm512_loadu_ps(&wrk2[i+128]);
                            _mm512_storeu_ps(&clogr[i+128],zmm25);
                            _mm512_storeu_ps(&clogi[i+128],zmm26);
                            zmm27 = _mm512_loadu_ps(&wrk1[i+144]);
                            zmm28 = xlogf(zmm27);
                            zmm29 = _mm512_loadu_ps(&wrk2[i+144]);
                            _mm512_storeu_ps(&clogr[i+144],zmm28);
                            _mm512_storeu_ps(&clogi[i+144],zmm29);
                            _mm_prefetch((const char*)&wrk1[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+192],_MM_HINT_T0);
                            zmm30 = _mm512_loadu_ps(&wrk1[i+160]);
                            zmm31 = xlogf(zmm30);
                            zmm0 = _mm512_loadu_ps(&wrk2[i+160]);
                            _mm512_storeu_ps(&clogr[i+160],zmm31);
                            _mm512_storeu_ps(&clogi[i+160],zmm0);
                            zmm1 = _mm512_loadu_ps(&wrk1[i+176]);
                            zmm2 = xlogf(zmm1);
                            zmm3 = _mm512_loadu_ps(&wrk2[i+176]);
                            _mm512_storeu_ps(&clogr[i+176],zmm2);
                            _mm512_storeu_ps(&clogi[i+176],zmm3);
                            _mm_prefetch((const char*)&wrk1[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+224],_MM_HINT_T0);
                            zmm4 = _mm512_loadu_ps(&wrk1[i+192]);
                            zmm5 = xlogf(zmm4);
                            zmm6 = _mm512_loadu_ps(&wrk2[i+192]);
                            _mm512_storeu_ps(&clogr[i+192],zmm5);
                            _mm512_storeu_ps(&clogi[i+192],zmm6);
                            zmm7 = _mm512_loadu_ps(&wrk1[i+208]);
                            zmm8 = xlogf(zmm7);
                            zmm9 = _mm512_loadu_ps(&wrk2[i+208]);
                            _mm512_storeu_ps(&clogr[i+208],zmm8);
                            _mm512_storeu_ps(&clogi[i+208],zmm9);
                            _mm_prefetch((const char*)&re[i+256],_MM_HINT_T0);
                            _mm_prefetch((const char*)&im[i+256],_MM_HINT_T0);
                            zmm10 = _mm512_loadu_ps(&wrk1[i+224]);
                            zmm11 = xlogf(zmm10);
                            zmm12 = _mm512_loadu_ps(&wrk2[i+224]);
                            _mm512_storeu_ps(&clogr[i+224],zmm11);
                            _mm512_storeu_ps(&clogi[i+224],zmm12);
                            zmm13 = _mm512_loadu_ps(&wrk1[i+240]);
                            zmm14 = xlogf(zmm13);
                            zmm15 = _mm512_loadu_ps(&wrk2[i+240]);
                            _mm512_storeu_ps(&clogr[i+240],zmm14);
                            _mm512_storeu_ps(&clogi[i+240],zmm15);
                       }

                        for(; (i+192) < n; i += 192) {
                             zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                            _mm512_storeu_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_loadu_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_storeu_ps(&clogr[i+80],zmm16);
                            _mm512_storeu_ps(&clogi[i+80],zmm17);
                            zmm18 = _mm512_loadu_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_loadu_ps(&wrk2[i+96]);
                            _mm512_storeu_ps(&clogr[i+96],zmm19);
                            _mm512_storeu_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_loadu_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_loadu_ps(&wrk2[i+112]);
                            _mm512_storeu_ps(&clogr[i+112],zmm22);
                            _mm512_storeu_ps(&clogi[i+112],zmm23);
                            zmm24 = _mm512_loadu_ps(&wrk1[i+128]);
                            zmm25 = xlogf(zmm24);
                            zmm26 = _mm512_loadu_ps(&wrk2[i+128]);
                            _mm512_storeu_ps(&clogr[i+128],zmm25);
                            _mm512_storeu_ps(&clogi[i+128],zmm26);
                            zmm27 = _mm512_loadu_ps(&wrk1[i+144]);
                            zmm28 = xlogf(zmm27);
                            zmm29 = _mm512_loadu_ps(&wrk2[i+144]);
                            _mm512_storeu_ps(&clogr[i+144],zmm28);
                            _mm512_storeu_ps(&clogi[i+144],zmm29);
                            zmm30 = _mm512_loadu_ps(&wrk1[i+160]);
                            zmm31 = xlogf(zmm30);
                            zmm0 = _mm512_loadu_ps(&wrk2[i+160]);
                            _mm512_storeu_ps(&clogr[i+160],zmm31);
                            _mm512_storeu_ps(&clogi[i+160],zmm0);
                            zmm1 = _mm512_loadu_ps(&wrk1[i+176]);
                            zmm2 = xlogf(zmm1);
                            zmm3 = _mm512_loadu_ps(&wrk2[i+176]);
                            _mm512_storeu_ps(&clogr[i+176],zmm2);
                            _mm512_storeu_ps(&clogi[i+176],zmm3);
                       }

                        for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                            _mm512_storeu_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_loadu_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_storeu_ps(&clogr[i+80],zmm16);
                            _mm512_storeu_ps(&clogi[i+80],zmm17);
                            zmm18 = _mm512_loadu_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_loadu_ps(&wrk2[i+96]);
                            _mm512_storeu_ps(&clogr[i+96],zmm19);
                            _mm512_storeu_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_loadu_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_loadu_ps(&wrk2[i+112]);
                            _mm512_storeu_ps(&clogr[i+112],zmm22);
                            _mm512_storeu_ps(&clogi[i+112],zmm23);
                       }

                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                       }

                        for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                       }

                        for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                       }

                         for(; (i+0) < n; i += 1) {
                            const float re   = wrk1[i];
                            const float x    = ceph_logf(re);
                            const float im   = wrk2[i];
                            clogr[i]         = x;
                            clogi[i]         = im; 
                       }
 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void clogv_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict __ATTR_ALIGN__(64) wrk1,
                                                   float * __restrict __ATTR_ALIGN__(64) wrk2,
                                                   float * __restrict __ATTR_ALIGN__(64) clogr,
                                                   float * __restrict __ATTR_ALIGN__(64) clogi,
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
                          cargv_zmm16r4_unroll_16x_a(re,im,wrk1,n);
                          cabsv_zmm16r4_unroll_16x_a(re,im,wrk2,n);
                          // wrk1 -- cabs
                          // wrk2 -- carg
                        for(; (i+255) < n; i += 256) {
                            _mm_prefetch((const char*)&wrk1[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            _mm_prefetch((const char*)&wrk1[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+64],_MM_HINT_T0);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            _mm_prefetch((const char*)&wrk1[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                            _mm512_store_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_load_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_load_ps(&wrk2[i+80]);
                            _mm512_store_ps(&clogr[i+80],zmm16);
                            _mm512_store_ps(&clogi[i+80],zmm17);
                            _mm_prefetch((const char*)&wrk1[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+128],_MM_HINT_T0);
                            zmm18 = _mm512_load_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_load_ps(&wrk2[i+96]);
                            _mm512_store_ps(&clogr[i+96],zmm19);
                            _mm512_store_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_load_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_load_ps(&wrk2[i+112]);
                            _mm512_store_ps(&clogr[i+112],zmm22);
                            _mm512_store_ps(&clogi[i+112],zmm23);
                            _mm_prefetch((const char*)&wrk1[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+160],_MM_HINT_T0);
                            zmm24 = _mm512_load_ps(&wrk1[i+128]);
                            zmm25 = xlogf(zmm24);
                            zmm26 = _mm512_load_ps(&wrk2[i+128]);
                            _mm512_store_ps(&clogr[i+128],zmm25);
                            _mm512_store_ps(&clogi[i+128],zmm26);
                            zmm27 = _mm512_load_ps(&wrk1[i+144]);
                            zmm28 = xlogf(zmm27);
                            zmm29 = _mm512_load_ps(&wrk2[i+144]);
                            _mm512_store_ps(&clogr[i+144],zmm28);
                            _mm512_store_ps(&clogi[i+144],zmm29);
                            _mm_prefetch((const char*)&wrk1[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+192],_MM_HINT_T0);
                            zmm30 = _mm512_load_ps(&wrk1[i+160]);
                            zmm31 = xlogf(zmm30);
                            zmm0 = _mm512_load_ps(&wrk2[i+160]);
                            _mm512_store_ps(&clogr[i+160],zmm31);
                            _mm512_store_ps(&clogi[i+160],zmm0);
                            zmm1 = _mm512_load_ps(&wrk1[i+176]);
                            zmm2 = xlogf(zmm1);
                            zmm3 = _mm512_load_ps(&wrk2[i+176]);
                            _mm512_store_ps(&clogr[i+176],zmm2);
                            _mm512_store_ps(&clogi[i+176],zmm3);
                            _mm_prefetch((const char*)&wrk1[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+224],_MM_HINT_T0);
                            zmm4 = _mm512_load_ps(&wrk1[i+192]);
                            zmm5 = xlogf(zmm4);
                            zmm6 = _mm512_load_ps(&wrk2[i+192]);
                            _mm512_store_ps(&clogr[i+192],zmm5);
                            _mm512_store_ps(&clogi[i+192],zmm6);
                            zmm7 = _mm512_load_ps(&wrk1[i+208]);
                            zmm8 = xlogf(zmm7);
                            zmm9 = _mm512_load_ps(&wrk2[i+208]);
                            _mm512_store_ps(&clogr[i+208],zmm8);
                            _mm512_store_ps(&clogi[i+208],zmm9);
                            _mm_prefetch((const char*)&re[i+256],_MM_HINT_T0);
                            _mm_prefetch((const char*)&im[i+256],_MM_HINT_T0);
                            zmm10 = _mm512_load_ps(&wrk1[i+224]);
                            zmm11 = xlogf(zmm10);
                            zmm12 = _mm512_load_ps(&wrk2[i+224]);
                            _mm512_store_ps(&clogr[i+224],zmm11);
                            _mm512_store_ps(&clogi[i+224],zmm12);
                            zmm13 = _mm512_load_ps(&wrk1[i+240]);
                            zmm14 = xlogf(zmm13);
                            zmm15 = _mm512_load_ps(&wrk2[i+240]);
                            _mm512_store_ps(&clogr[i+240],zmm14);
                            _mm512_store_ps(&clogi[i+240],zmm15);
                       }

                        for(; (i+192) < n; i += 192) {
                             zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                            _mm512_store_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_load_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_load_ps(&wrk2[i+80]);
                            _mm512_store_ps(&clogr[i+80],zmm16);
                            _mm512_store_ps(&clogi[i+80],zmm17);
                            zmm18 = _mm512_load_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_load_ps(&wrk2[i+96]);
                            _mm512_store_ps(&clogr[i+96],zmm19);
                            _mm512_store_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_load_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_load_ps(&wrk2[i+112]);
                            _mm512_store_ps(&clogr[i+112],zmm22);
                            _mm512_store_ps(&clogi[i+112],zmm23);
                            zmm24 = _mm512_load_ps(&wrk1[i+128]);
                            zmm25 = xlogf(zmm24);
                            zmm26 = _mm512_load_ps(&wrk2[i+128]);
                            _mm512_store_ps(&clogr[i+128],zmm25);
                            _mm512_store_ps(&clogi[i+128],zmm26);
                            zmm27 = _mm512_load_ps(&wrk1[i+144]);
                            zmm28 = xlogf(zmm27);
                            zmm29 = _mm512_load_ps(&wrk2[i+144]);
                            _mm512_store_ps(&clogr[i+144],zmm28);
                            _mm512_store_ps(&clogi[i+144],zmm29);
                            zmm30 = _mm512_load_ps(&wrk1[i+160]);
                            zmm31 = xlogf(zmm30);
                            zmm0 = _mm512_load_ps(&wrk2[i+160]);
                            _mm512_store_ps(&clogr[i+160],zmm31);
                            _mm512_store_ps(&clogi[i+160],zmm0);
                            zmm1 = _mm512_load_ps(&wrk1[i+176]);
                            zmm2 = xlogf(zmm1);
                            zmm3 = _mm512_load_ps(&wrk2[i+176]);
                            _mm512_store_ps(&clogr[i+176],zmm2);
                            _mm512_store_ps(&clogi[i+176],zmm3);
                       }

                        for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                            _mm512_store_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_load_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_store_ps(&clogr[i+80],zmm16);
                            _mm512_store_ps(&clogi[i+80],zmm17);
                            zmm18 = _mm512_load_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_load_ps(&wrk2[i+96]);
                            _mm512_store_ps(&clogr[i+96],zmm19);
                            _mm512_store_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_load_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_load_ps(&wrk2[i+112]);
                            _mm512_store_ps(&clogr[i+112],zmm22);
                            _mm512_store_ps(&clogi[i+112],zmm23);
                       }

                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                       }

                        for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                       }

                        for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                       }

                         for(; (i+0) < n; i += 1) {
                            const float re   = wrk1[i];
                            const float x    = ceph_logf(re);
                            const float im   = wrk2[i];
                            clogr[i]         = x;
                            clogi[i]         = im; 
                       }
 
                 }


                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void clogv_zmm16r4_unroll_10x_u(const float * __restrict re,
                                                   const float * __restrict im,
                                                   float * __restrict wrk1,
                                                   float * __restrict wrk2,
                                                   float * __restrict clogr,
                                                   float * __restrict clogi,
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
                          cargv_zmm16r4_unroll_10x_u(re,im,wrk1,n);
                          cabsv_zmm16r4_unroll_10x_u(re,im,wrk2,n);
                          // wrk1 -- cabs
                          // wrk2 -- carg
                        for(; (i+159) < n; i += 160) {
                            _mm_prefetch((const char*)&wrk1[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            _mm_prefetch((const char*)&wrk1[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+64],_MM_HINT_T0);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            _mm_prefetch((const char*)&wrk1[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                            _mm512_storeu_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_loadu_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_storeu_ps(&clogr[i+80],zmm16);
                            _mm512_storeu_ps(&clogi[i+80],zmm17);
                            _mm_prefetch((const char*)&wrk1[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+128],_MM_HINT_T0);
                            zmm18 = _mm512_loadu_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_loadu_ps(&wrk2[i+96]);
                            _mm512_storeu_ps(&clogr[i+96],zmm19);
                            _mm512_storeu_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_loadu_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_loadu_ps(&wrk2[i+112]);
                            _mm512_storeu_ps(&clogr[i+112],zmm22);
                            _mm512_storeu_ps(&clogi[i+112],zmm23);
                            _mm_prefetch((const char*)&wrk1[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+160],_MM_HINT_T0);
                            zmm24 = _mm512_loadu_ps(&wrk1[i+128]);
                            zmm25 = xlogf(zmm24);
                            zmm26 = _mm512_loadu_ps(&wrk2[i+128]);
                            _mm512_storeu_ps(&clogr[i+128],zmm25);
                            _mm512_storeu_ps(&clogi[i+128],zmm26);
                            zmm27 = _mm512_loadu_ps(&wrk1[i+144]);
                            zmm28 = xlogf(zmm27);
                            zmm29 = _mm512_loadu_ps(&wrk2[i+144]);
                            _mm512_storeu_ps(&clogr[i+144],zmm28);
                            _mm512_storeu_ps(&clogi[i+144],zmm29);
                          
                       }

                       

                        for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                            _mm512_storeu_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_loadu_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_storeu_ps(&clogr[i+80],zmm16);
                            _mm512_storeu_ps(&clogi[i+80],zmm17);
                            zmm18 = _mm512_loadu_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_loadu_ps(&wrk2[i+96]);
                            _mm512_storeu_ps(&clogr[i+96],zmm19);
                            _mm512_storeu_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_loadu_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_loadu_ps(&wrk2[i+112]);
                            _mm512_storeu_ps(&clogr[i+112],zmm22);
                            _mm512_storeu_ps(&clogi[i+112],zmm23);
                       }

                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                       }

                        for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                       }

                        for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                       }

                         for(; (i+0) < n; i += 1) {
                            const float re   = wrk1[i];
                            const float x    = ceph_logf(re);
                            const float im   = wrk2[i];
                            clogr[i]         = x;
                            clogi[i]         = im; 
                       }
 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void clogv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict __ATTR_ALIGN__(64) wrk1,
                                                   float * __restrict __ATTR_ALIGN__(64) wrk2,
                                                   float * __restrict __ATTR_ALIGN__(64) clogr,
                                                   float * __restrict __ATTR_ALIGN__(64) clogi,
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
                          cargv_zmm16r4_unroll_10x_a(re,im,wrk1,n);
                          cabsv_zmm16r4_unroll_10x_a(re,im,wrk2,n);
                          // wrk1 -- cabs
                          // wrk2 -- carg
                        for(; (i+159) < n; i += 160) {
                            _mm_prefetch((const char*)&wrk1[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            _mm_prefetch((const char*)&wrk1[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+64],_MM_HINT_T0);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            _mm_prefetch((const char*)&wrk1[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                            _mm512_store_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_load_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_load_ps(&wrk2[i+80]);
                            _mm512_store_ps(&clogr[i+80],zmm16);
                            _mm512_store_ps(&clogi[i+80],zmm17);
                            _mm_prefetch((const char*)&wrk1[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+128],_MM_HINT_T0);
                            zmm18 = _mm512_load_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_load_ps(&wrk2[i+96]);
                            _mm512_store_ps(&clogr[i+96],zmm19);
                            _mm512_store_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_load_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_load_ps(&wrk2[i+112]);
                            _mm512_store_ps(&clogr[i+112],zmm22);
                            _mm512_store_ps(&clogi[i+112],zmm23);
                            _mm_prefetch((const char*)&wrk1[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+160],_MM_HINT_T0);
                            zmm24 = _mm512_load_ps(&wrk1[i+128]);
                            zmm25 = xlogf(zmm24);
                            zmm26 = _mm512_load_ps(&wrk2[i+128]);
                            _mm512_store_ps(&clogr[i+128],zmm25);
                            _mm512_store_ps(&clogi[i+128],zmm26);
                            zmm27 = _mm512_load_ps(&wrk1[i+144]);
                            zmm28 = xlogf(zmm27);
                            zmm29 = _mm512_load_ps(&wrk2[i+144]);
                            _mm512_store_ps(&clogr[i+144],zmm28);
                            _mm512_store_ps(&clogi[i+144],zmm29);
                           
                       }

                       
                        for(; (i+127) < n; i += 128) {
                             zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                            _mm512_store_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_load_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_store_ps(&clogr[i+80],zmm16);
                            _mm512_store_ps(&clogi[i+80],zmm17);
                            zmm18 = _mm512_load_ps(&wrk1[i+96]);
                            zmm19 = xlogf(zmm18);
                            zmm20 = _mm512_load_ps(&wrk2[i+96]);
                            _mm512_store_ps(&clogr[i+96],zmm19);
                            _mm512_store_ps(&clogi[i+96],zmm20);
                            zmm21 = _mm512_load_ps(&wrk1[i+112]);
                            zmm22 = xlogf(zmm21);
                            zmm23 = _mm512_load_ps(&wrk2[i+112]);
                            _mm512_store_ps(&clogr[i+112],zmm22);
                            _mm512_store_ps(&clogi[i+112],zmm23);
                       }

                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                       }

                        for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                       }

                        for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                       }

                         for(; (i+0) < n; i += 1) {
                            const float re   = wrk1[i];
                            const float x    = ceph_logf(re);
                            const float im   = wrk2[i];
                            clogr[i]         = x;
                            clogi[i]         = im; 
                       }
 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void clogv_zmm16r4_unroll_6x_u(const float * __restrict re,
                                                   const float * __restrict im,
                                                   float * __restrict wrk1,
                                                   float * __restrict wrk2,
                                                   float * __restrict clogr,
                                                   float * __restrict clogi,
                                                   const int32_t n) {

                        if(__builtin_expect(0==n,0)) { return;}
                          register vfloat zmm0,zmm1,zmm2,zmm3;
                          register vfloat zmm4,zmm5,zmm6,zmm7;
                          register vfloat zmm8,zmm9,zmm10,zmm11;
                          register vfloat zmm12,zmm13,zmm14,zmm15;
                          register vfloat zmm16,zmm17;
                          int32_t i;
                          cargv_zmm16r4_unroll_6x_u(re,im,wrk1,n);
                          cabsv_zmm16r4_unroll_6x_u(re,im,wrk2,n);
                          // wrk1 -- cabs
                          // wrk2 -- carg
                        for(; (i+95) < n; i += 96) {
                            _mm_prefetch((const char*)&wrk1[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            _mm_prefetch((const char*)&wrk1[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+64],_MM_HINT_T0);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            _mm_prefetch((const char*)&wrk1[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                            _mm512_storeu_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_loadu_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_loadu_ps(&wrk2[i+80]);
                            _mm512_storeu_ps(&clogr[i+80],zmm16);
                            _mm512_storeu_ps(&clogi[i+80],zmm17);
                                                     
                       }

                       for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_loadu_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_loadu_ps(&wrk2[i+64]);
                            _mm512_storeu_ps(&clogr[i+64],zmm13);
                       }

                        for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_loadu_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_loadu_ps(&wrk2[i+32]);
                            _mm512_storeu_ps(&clogr[i+32],zmm7);
                            _mm512_storeu_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_loadu_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_loadu_ps(&wrk2[i+48]);
                            _mm512_storeu_ps(&clogr[i+48],zmm10);
                            _mm512_storeu_ps(&clogi[i+48],zmm11);
                       }

                        for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_loadu_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_storeu_ps(&clogr[i+16],zmm4);
                            _mm512_storeu_ps(&clogi[i+16],zmm5);
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_loadu_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_loadu_ps(&wrk2[i+0]);
                            _mm512_storeu_ps(&clogr[i+0],zmm1);
                            _mm512_storeu_ps(&clogi[i+0],zmm2);
                       }

                         for(; (i+0) < n; i += 1) {
                            const float re   = wrk1[i];
                            const float x    = ceph_logf(re);
                            const float im   = wrk2[i];
                            clogr[i]         = x;
                            clogi[i]         = im; 
                       }
 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void clogv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict __ATTR_ALIGN__(64) wrk1,
                                                   float * __restrict __ATTR_ALIGN__(64) wrk2,
                                                   float * __restrict __ATTR_ALIGN__(64) clogr,
                                                   float * __restrict __ATTR_ALIGN__(64) clogi,
                                                   const int32_t n) {

                        if(__builtin_expect(0==n,0)) { return;}
                          register vfloat zmm0,zmm1,zmm2,zmm3;
                          register vfloat zmm4,zmm5,zmm6,zmm7;
                          register vfloat zmm8,zmm9,zmm10,zmm11;
                          register vfloat zmm12,zmm13,zmm14,zmm15;
                          register vfloat zmm16,zmm17;
                          int32_t i;
                          cargv_zmm16r4_unroll_6x_a(re,im,wrk1,n);
                          cabsv_zmm16r4_unroll_6x_a(re,im,wrk2,n);
                          // wrk1 -- cabs
                          // wrk2 -- carg
                        for(; (i+95) < n; i += 96) {
                            _mm_prefetch((const char*)&wrk1[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+32],_MM_HINT_T0);
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_loadu_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            _mm_prefetch((const char*)&wrk1[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+64],_MM_HINT_T0);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            _mm_prefetch((const char*)&wrk1[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&wrk2[i+96],_MM_HINT_T0);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                            _mm512_store_ps(&clogi[i+64],zmm14);
                            zmm15 = _mm512_load_ps(&wrk1[i+80]);
                            zmm16 = xlogf(zmm15);
                            zmm17 = _mm512_load_ps(&wrk2[i+80]);
                            _mm512_store_ps(&clogr[i+80],zmm16);
                            _mm512_store_ps(&clogi[i+80],zmm17);
                                                      
                       }

                        for(; (i+79) < n; i += 80) {
                             zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                            zmm12 = _mm512_load_ps(&wrk1[i+64]);
                            zmm13 = xlogf(zmm12);
                            zmm14 = _mm512_load_ps(&wrk2[i+64]);
                            _mm512_store_ps(&clogr[i+64],zmm13);
                       }

                        for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                            zmm6 = _mm512_load_ps(&wrk1[i+32]);
                            zmm7 = xlogf(zmm6);
                            zmm8 = _mm512_load_ps(&wrk2[i+32]);
                            _mm512_store_ps(&clogr[i+32],zmm7);
                            _mm512_store_ps(&clogi[i+32],zmm8);
                            zmm9 = _mm512_load_ps(&wrk1[i+48]);
                            zmm10= xlogf(zmm9);
                            zmm11= _mm512_load_ps(&wrk2[i+48]);
                            _mm512_store_ps(&clogr[i+48],zmm10);
                            _mm512_store_ps(&clogi[i+48],zmm11);
                       }

                        for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                            zmm3 = _mm512_load_ps(&wrk1[i+16]);
                            zmm4 = xlogf(zmm3);
                            zmm5 = _mm512_load_ps(&wrk2[i+16]);
                            _mm512_store_ps(&clogr[i+16],zmm4);
                            _mm512_store_ps(&clogi[i+16],zmm5);
                       }

                         for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_load_ps(&wrk1[i+0]);
                            zmm1 = xlogf(zmm0);
                            zmm2 = _mm512_load_ps(&wrk2[i+0]);
                            _mm512_store_ps(&clogr[i+0],zmm1);
                            _mm512_store_ps(&clogi[i+0],zmm2);
                       }

                         for(; (i+0) < n; i += 1) {
                            const float re   = wrk1[i];
                            const float x    = ceph_logf(re);
                            const float im   = wrk2[i];
                            clogr[i]         = x;
                            clogi[i]         = im; 
                       }
 
                 }



       } // math



} // gms









#endif /*__GMS_CLOG_VEC_ZMM16R4_HPP__*/
