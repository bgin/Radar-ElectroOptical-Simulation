

#ifndef __GMS_RCS_CYLINDER_ZMM16R4_OMP_HPP__
#define __GMS_RCS_CYLINDER_ZMM16R4_OMP_HPP__

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

    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_MAJOR = 1U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_MINOR = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_MICRO = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_FULLVER =
      1000U*GMS_RCS_CYLINDER_ZMM16R4_OMP_MAJOR+
      100U*GMS_RCS_CYLINDER_ZMM16R4_OMP_MINOR+
      10U*GMS_RCS_CYLINDER_ZMM16R4_OMP_MICRO;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_CREATION_DATE = "20-01-2023 16:36 PM +00200 (FRI 20 JAN 2023 GMT+2)";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_DESCRIPTION   = "AVX512 optimized Cylinder Radar Cross Section (analytic) functionality OpenMP accelerated.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_rcs_cylinder_zmm16r4.hpp"


namespace gms {


           namespace radiolocation {
           
                        
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f419_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 16;
                         __m512 a1,a2,a3,a4,a5,a6,a7;
                         __m512 a8,a9,a10,a11,a12,a13,a14,a15,a16;
                         __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8;
                         __m512 k0a9,k0a10,k0a11,k0a12,k0a13,k0a14,k0a15,k0a16;
                         __m512 rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8;
                         __m512 rcs9,rcs10,rcs11,rcs12,rcs13,rcs14,rcs15,rcs16;
                         int32_t j,m,m1;
                         
                         m = n%16;
                         if(m!=0) {
                            for(j = 0; j != m; ++j) {
                                a   = pa[j];
                                k0a = pk0a[j];
                                rcs = rcs_f419_zmm16r4(a,k0a);
                                prcs[j] = rcs;
                            }
                            if(n<16) {return;}
                         }                    
                         
                         m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                              \
        firstprivate(m1,PF_DIST) private(j,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16)    \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8)                     \
                                 private(k0a9,k0a10,k0a11,k0a12,k0a13,k0a14,k0a15,k0a16)              \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8)                     \
                                 private(rcs9,rcs10,rcs11,rcs12,rcs13,rcs14,rcs15,rcs16)              \
                                 shared(n,pa,pk0a,prcs)
                         for(j = m1; j != n; j += 16) {
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif             
                              a1   = pa[j+0];
                              k0a1 = pk0a[j+0];
                              rcs = rcs_f419_zmm16r4(a1,k0a1);
                              prcs[j+0] = rcs1;  
                              a2   = pa[j+1];
                              k0a2 = pk0a[j+1];
                              rcs2 = rcs_f419_zmm16r4(a2,k0a2);
                              prcs[j+1] = rcs2;
                              a3   = pa[j+2];
                              k0a3 = pk0a[j+2];
                              rcs3 = rcs_f419_zmm16r4(a3,k0a3);
                              prcs[j+2] = rcs3; 
                              a4   = pa[j+3];
                              k0a4 = pk0a[j+3];
                              rcs4 = rcs_f419_zmm16r4(a4,k0a4);
                              prcs[j+3] = rcs4;   
                              a5   = pa[j+4];
                              k0a55 = pk0a[j+4];
                              rcs5 = rcs_f419_zmm16r4(a5,k0a5);
                              prcs[j+4] = rcs5;
                              a6   = pa[j+5];
                              k0a6 = pk0a[j+5];
                              rcs6 = rcs_f419_zmm16r4(a6,k0a6);
                              prcs[j+5] = rcs6; 
                              a7   = pa[j+6];
                              k0a7 = pk0a[j+6];
                              rcs7 = rcs_f419_zmm16r4(a7,k0a7);
                              prcs[j+6] = rcs7;   
                              a8   = pa[j+7];
                              k0a8 = pk0a[j+7];
                              rcs8 = rcs_f419_zmm16r4(a8,k0a8);
                              prcs[j+7] = rcs8;
                              a9   = pa[j+8];
                              k0a9 = pk0a[j+8];
                              rcs9 = rcs_f419_zmm16r4(a9,k0a9);
                              prcs[j+8] = rcs9;  
                              a10   = pa[j+9];
                              k0a10 = pk0a[j+9];
                              rcs10 = rcs_f419_zmm16r4(a10,k0a10);
                              prcs[j+9] = rcs10;
                              a11   = pa[j+10];
                              k0a11 = pk0a[j+10];
                              rcs11 = rcs_f419_zmm16r4(a11,k0a11);
                              prcs[j+10] = rcs11; 
                              a12   = pa[j+11];
                              k0a12 = pk0a[j+11];
                              rcs12 = rcs_f419_zmm16r4(a12,k0a12);
                              prcs[j+11] = rcs12;
                              a13   = pa[j+12];
                              k0a13 = pk0a[j+12];
                              rcs13 = rcs_f419_zmm16r4(a13,k0a13);
                              prcs[j+12] = rcs13;
                              a14   = pa[j+13];
                              k0a14 = pk0a[j+13];
                              rcs14 = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+13] = rcs14;
                              a15   = pa[j+14];
                              k0a15 = pk0a[j+14];
                              rcs15 = rcs_f419_zmm16r4(a15,k0a15);
                              prcs[j+14] = rcs15;
                              a16   = pa[j+15];
                              k0a16 = pk0a[j+15];
                              rcs16 = rcs_f419_zmm16r4(a16,k0a16);
                              prcs[j+15] = rcs16;
                         }            
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f419_zmm16r4_unroll10x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 10;
                         __m512 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10;
                         __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8,k0a9,k0a10;
                         __m512 rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8,rcs9,rcs10;
                         int32_t j,m,m1;
                         
                         m = n%10;
                         if(m!=0) {
                            for(j = 0; j != m; ++j) {
                                a   = pa[j];
                                k0a = pk0a[j];
                                rcs = rcs_f419_zmm16r4(a,k0a);
                                prcs[j] = rcs;
                            }
                            if(n<10) {return;}
                         }                    
                         
                         m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                             \
        firstprivate(m1,PF_DIST) private(j,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)                           \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8,k0a9,k0a10)         \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8,rcs8,rcs9,rcs10)    \
                                 shared(n,pa,pk0a,prcs)
                         for(j = m1; j != n; j += 10) {
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif             
                              a1   = pa[j+0];
                              k0a1 = pk0a[j+0];
                              rcs = rcs_f419_zmm16r4(a1,k0a1);
                              prcs[j+0] = rcs1;  
                              a2   = pa[j+1];
                              k0a2 = pk0a[j+1];
                              rcs2 = rcs_f419_zmm16r4(a2,k0a2);
                              prcs[j+1] = rcs2;
                              a3   = pa[j+2];
                              k0a3 = pk0a[j+2];
                              rcs3 = rcs_f419_zmm16r4(a3,k0a3);
                              prcs[j+2] = rcs3; 
                              a4   = pa[j+3];
                              k0a4 = pk0a[j+3];
                              rcs4 = rcs_f419_zmm16r4(a4,k0a4);
                              prcs[j+3] = rcs4;   
                              a5   = pa[j+4];
                              k0a55 = pk0a[j+4];
                              rcs5 = rcs_f419_zmm16r4(a5,k0a5);
                              prcs[j+4] = rcs5;
                              a6   = pa[j+5];
                              k0a6 = pk0a[j+5];
                              rcs6 = rcs_f419_zmm16r4(a6,k0a6);
                              prcs[j+5] = rcs6; 
                              a7   = pa[j+6];
                              k0a7 = pk0a[j+6];
                              rcs7 = rcs_f419_zmm16r4(a7,k0a7);
                              prcs[j+6] = rcs7;   
                              a8   = pa[j+7];
                              k0a8 = pk0a[j+7];
                              rcs8 = rcs_f419_zmm16r4(a8,k0a8);
                              prcs[j+7] = rcs8;
                              a9   = pa[j+8];
                              k0a9 = pk0a[j+8];
                              rcs9 = rcs_f419_zmm16r4(a9,k0a9);
                              prcs[j+8] = rcs9;  
                              a10   = pa[j+9];
                              k0a10 = pk0a[j+9];
                              rcs10 = rcs_f419_zmm16r4(a10,k0a10);
                              prcs[j+9] = rcs10;
                             
                         }            
              }
              
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f419_zmm16r4_unroll6x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 10;
                         __m512 a1,a2,a3,a4,a5,a6;
                         __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6;
                         __m512 rcs1,rcs2,rcs3,rcs4,rcs5,rcs6;
                         int32_t j,m,m1;
                         
                         m = n%6;
                         if(m!=0) {
                            for(j = 0; j != m; ++j) {
                                a   = pa[j];
                                k0a = pk0a[j];
                                rcs = rcs_f419_zmm16r4(a,k0a);
                                prcs[j] = rcs;
                            }
                            if(n<6) {return;}
                         }                    
                         
                         m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                   \
        firstprivate(m1,PF_DIST) private(j,a1,a2,a3,a4,a5,a6)              \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6)    \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6)    \
                                 shared(n,pa,pk0a,prcs)
                         for(j = m1; j != n; j += 6) {
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif             
                              a1   = pa[j+0];
                              k0a1 = pk0a[j+0];
                              rcs = rcs_f419_zmm16r4(a1,k0a1);
                              prcs[j+0] = rcs1;  
                              a2   = pa[j+1];
                              k0a2 = pk0a[j+1];
                              rcs2 = rcs_f419_zmm16r4(a2,k0a2);
                              prcs[j+1] = rcs2;
                              a3   = pa[j+2];
                              k0a3 = pk0a[j+2];
                              rcs3 = rcs_f419_zmm16r4(a3,k0a3);
                              prcs[j+2] = rcs3; 
                              a4   = pa[j+3];
                              k0a4 = pk0a[j+3];
                              rcs4 = rcs_f419_zmm16r4(a4,k0a4);
                              prcs[j+3] = rcs4;   
                              a5   = pa[j+4];
                              k0a55 = pk0a[j+4];
                              rcs5 = rcs_f419_zmm16r4(a5,k0a5);
                              prcs[j+4] = rcs5;
                              a6   = pa[j+5];
                              k0a6 = pk0a[j+5];
                              rcs6 = rcs_f419_zmm16r4(a6,k0a6);
                              prcs[j+5] = rcs6; 
                            
                             
                         }            
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f4120_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 16;
                         register __m512 a1,a2,a3,a4,a5,a6,a7,a8;
                         register __m512 a9,a10,a11,a12,a13,a14,a15,a16;
                         register __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8;
                         register __m512 k0a9,k0a10,k0a11,k0a12,k0a13,k0a14,k0a15,k0a16;
                         register __m512 rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8;
                         register __m512 rcs9,rcs10,rcs11,rcs12,rcs13,rcs14,rcs15,rcs16;
                         int32_t j,m,m1;
                         
                         m = n%16;
                         if(m!=0) {
                            for(j = 0; j != m; ++j) {
                                a   = pa[j];
                                k0a = pk0a[j];
                                rcs = rcs_f4120_zmm16r4(a,k0a);
                                prcs[j] = rcs;
                            }
                            if(n<16) {return;}
                         }                    
                         
                         m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                              \
        firstprivate(m1,PF_DIST) private(j,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16)    \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8)                     \
                                 private(k0a9,k0a10,k0a11,k0a12,k0a13,k0a14,k0a15,k0a16)              \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8)                     \
                                 private(rcs9,rcs10,rcs11,rcs12,rcs13,rcs14,rcs15,rcs16)              \
                                 shared(n,pa,pk0a,prcs)
                         for(j = m1; j != n; j += 16) {[
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif             
                              a   = pa[j+0];
                              k0a = pk0a[j+0];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+0] = rcs;  
                              a   = pa[j+1];
                              k0a = pk0a[j+1];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+1] = rcs;
                              a   = pa[j+2];
                              k0a = pk0a[j+2];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+2] = rcs; 
                              a   = pa[j+3];
                              k0a = pk0a[j+3];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+3] = rcs;   
                              a   = pa[j+4];
                              k0a = pk0a[j+4];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+4] = rcs;
                              a   = pa[j+5];
                              k0a = pk0a[j+5];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+5] = rcs; 
                              a   = pa[j+6];
                              k0a = pk0a[j+6];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+6] = rcs;   
                              a   = pa[j+7];
                              k0a = pk0a[j+7];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+7] = rcs;
                              a   = pa[j+8];
                              k0a = pk0a[j+8];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+8] = rcs;  
                              a   = pa[j+9];
                              k0a = pk0a[j+9];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+9] = rcs;
                              a   = pa[j+10];
                              k0a = pk0a[j+10];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+10] = rcs; 
                              a   = pa[j+11];
                              k0a = pk0a[j+11];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+11] = rcs;
                              a   = pa[j+12];
                              k0a = pk0a[j+12];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+12] = rcs;
                              a   = pa[j+13];
                              k0a = pk0a[j+13];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+13] = rcs;
                              a   = pa[j+14];
                              k0a = pk0a[j+14];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+14] = rcs;
                              a   = pa[j+15];
                              k0a = pk0a[j+15];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+15] = rcs;
                         }            
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f4120_zmm16r4_unroll10x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 10;
                         register __m512 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10;
                         register __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8,k0a9,k0a10;
                         register __m512 rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8,rcs9,rcs10;
                         int32_t j,m,m1;
                         
                         m = n%10;
                         if(m!=0) {
                            for(j = 0; j != m; ++j) {
                                a   = pa[j];
                                k0a = pk0a[j];
                                rcs = rcs_f4120_zmm16r4(a,k0a);
                                prcs[j] = rcs;
                            }
                            if(n<10) {return;}
                         }                    
                         
                         m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                              \
        firstprivate(m1,PF_DIST) private(j,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)                            \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8,k0a9,k0a10)          \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8,rcs9,rcs10)          \
                                 shared(n,pa,pk0a,prcs)
                         for(j = m1; j != n; j += 10) {
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif             
                              a   = pa[j+0];
                              k0a = pk0a[j+0];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+0] = rcs;  
                              a   = pa[j+1];
                              k0a = pk0a[j+1];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+1] = rcs;
                              a   = pa[j+2];
                              k0a = pk0a[j+2];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+2] = rcs; 
                              a   = pa[j+3];
                              k0a = pk0a[j+3];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+3] = rcs;   
                              a   = pa[j+4];
                              k0a = pk0a[j+4];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+4] = rcs;
                              a   = pa[j+5];
                              k0a = pk0a[j+5];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+5] = rcs; 
                              a   = pa[j+6];
                              k0a = pk0a[j+6];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+6] = rcs;   
                              a   = pa[j+7];
                              k0a = pk0a[j+7];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+7] = rcs;
                              a   = pa[j+8];
                              k0a = pk0a[j+8];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+8] = rcs;  
                              a   = pa[j+9];
                              k0a = pk0a[j+9];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+9] = rcs;
                            
                         }            
              }
              
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f4120_zmm16r4_unroll6x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 6;
                         register __m512 a1,a2,a3,a4,a5,a6;
                         register __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6;
                         register __m512 rcs1,rcs2,rcs3,rcs4,rcs5,rcs6;
                         int32_t j,m,m1;
                         
                         m = n%6;
                         if(m!=0) {
                            for(j = 0; j != m; ++j) {
                                a   = pa[j];
                                k0a = pk0a[j];
                                rcs = rcs_f4120_zmm16r4(a,k0a);
                                prcs[j] = rcs;
                            }
                            if(n<6) {return;}
                         }                    
                         
                         m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                         \
        firstprivate(m1,PF_DIST) private(j,a1,a2,a3,a4,a5,a6)                    \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6)          \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6)          \
                                 shared(n,pa,pk0a,prcs)
                         for(j = m1; j != n; j += 6) {
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif             
                              a   = pa[j+0];
                              k0a = pk0a[j+0];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+0] = rcs;  
                              a   = pa[j+1];
                              k0a = pk0a[j+1];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+1] = rcs;
                              a   = pa[j+2];
                              k0a = pk0a[j+2];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+2] = rcs; 
                              a   = pa[j+3];
                              k0a = pk0a[j+3];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+3] = rcs;   
                              a   = pa[j+4];
                              k0a = pk0a[j+4];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+4] = rcs;
                              a   = pa[j+5];
                              k0a = pk0a[j+5];
                              rcs = rcs_f4120_zmm16r4(a,k0a);
                              prcs[j+5] = rcs; 
                                                        
                         }            
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void rcs_f4122_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pphi,
	                                            const __m512 * __restrict __ATTR_ALIGN__(64) pa,
	                                            const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
	                                            __m512 * __restrict __ATTR_ALIGN__(64) prcs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                            
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 16;
	                register __m512 phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8;
	                register __m512 phi9,phi10,phi11,phi12,phi13,phi14,phi15,phi16;
	                register __m512 a1,a2,a3,a4,a5,a6,a7,a8;
	                register __m512 a9,a10,a11,a12,a13,a14,a15,a16;
	                register __m512 k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8;
	                register __m512 k0a9,k0a10,k0a11,k0a12,k0a13,k0a14,k0a15,k0a16;
	                register __m512 rcs1,rcs2,rcs3,rcs4,rcs,rcs6,rcs7,rcs8;
	                register __m512 rcs9,rcs10,rcs11,rcs12,rcs13,rcs14,rcs15,rcs16;
	                int32_t j,m,m1;
	                
	                m = n%16;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       phi = phi[j];
	                       a   = pa[j];
	                       k0a = pk0a[j];
	                       rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                       prcs[j] = rcs;
	                   }
	                   if(n<16) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                              \
        firstprivate(m1,PF_DIST) private(j,phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8)                   \
                                 private(phi9,phi10,phi11,phi12,phi13,phi14,phi15,phi16)              \
                                 private(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16)      \
                                 private(k0a1,k0a2,k0a3,k0a4,k0a5,k0a6,k0a7,k0a8)                     \
                                 private(k0a9,k0a10,k0a11,k0a12,k0a13,k0a14,k0a15,k0a16)              \
                                 private(rcs1,rcs2,rcs3,rcs4,rcs5,rcs6,rcs7,rcs8)                     \
                                 private(rcs9,rcs10,rcs11,rcs12,rcs13,rcs14,rcs15,rcs16)              \
                                 shared(n,pa,pk0a,prcs)	                
	                for(j = m1; j != n; j += 16) {
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                              _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T0);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T0);
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2
                              _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T1);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T1);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 3
                              _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_T2);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_T2);
#elif  (__RCS_CYLINDER_PF_CACHE_HINT__) == 4
                              _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pa[j+PF_DIST],_MM_HINT_NTA);
                              _mm_prefetch((char*)&pk0a[j+PF_DIST],_MM_HINT_NTA);
#endif  	              
                              phi = phi[j+0];
	                      a   = pa[j+0];
	                      k0a = pk0a[j+0];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+0] = rcs; 
	                      phi = phi[j+1];
	                      a   = pa[j+1];
	                      k0a = pk0a[j+1];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+1] = rcs;  
	                      phi = phi[j+2];
	                      a   = pa[j+2];
	                      k0a = pk0a[j+2];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+2] = rcs; 
	                      phi = phi[j+3];
	                      a   = pa[j+3];
	                      k0a = pk0a[j+3];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+3] = rcs;
	                      phi = phi[j+4];
	                      a   = pa[j+4];
	                      k0a = pk0a[j+4];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+4] = rcs;
	                      phi = phi[j+5];
	                      a   = pa[j+5];
	                      k0a = pk0a[j+5];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+5] = rcs;
	                      phi = phi[j+6];
	                      a   = pa[j+6];
	                      k0a = pk0a[j+6];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+6] = rcs;
	                      phi = phi[j+7];
	                      a   = pa[j+7];
	                      k0a = pk0a[j+7];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+7] = rcs;
	                      phi = phi[j+8];
	                      a   = pa[j+8];
	                      k0a = pk0a[j+8];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+8] = rcs;
	                      phi = phi[j+9];
	                      a   = pa[j+9];
	                      k0a = pk0a[j+9];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+9] = rcs;
	                      phi = phi[j+10];
	                      a   = pa[j+10];
	                      k0a = pk0a[j+10];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+10] = rcs;
	                      phi = phi[j+11];
	                      a   = pa[j+11];
	                      k0a = pk0a[j+11];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+11] = rcs; 
	                      phi = phi[j+12];
	                      a   = pa[j+12];
	                      k0a = pk0a[j+12];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+12] = rcs;  
	                      phi = phi[j+13];
	                      a   = pa[j+13];
	                      k0a = pk0a[j+13];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+13] = rcs; 
	                      phi = phi[j+14];
	                      a   = pa[j+14];
	                      k0a = pk0a[j+14];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+14] = rcs;
	                      phi = phi[j+15];
	                      a   = pa[j+15];
	                      k0a = pk0a[j+15];
	                      rcs = rcs_f4122_zmm16r4(phi,a,k0a);
	                      prcs[j+15] = rcs;         
	                }              
	       } 
	       
              
              
              
              
              
           
         } // radiolocation

} // gms










#endif /*__GMS_RCS_CYLINDER_ZMM16R4_OMP_HPP__*/
