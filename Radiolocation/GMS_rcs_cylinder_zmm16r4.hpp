
#ifndef __GMS_RCS_CYLINDER_ZMM16R4_HPP__
#define __GMS_RCS_CYLINDER_ZMM16R4_HPP__ 200120231636


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

    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_FULLVER =
      1000U*GMS_RCS_CYLINDER_ZMM16R4_MAJOR+
      100U*GMS_RCS_CYLINDER_ZMM16R4_MINOR+
      10U*GMS_RCS_CYLINDER_ZMM16R4_MICRO;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_CREATION_DATE = "20-01-2023 16:36 PM +00200 (FRI 20 JAN 2023 GMT+2)";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_DESCRIPTION   = "AVX512 optimized Cylinder Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"


#ifndef __RCS_CYLINDER_PF_CACHE_HINT__
#define __RCS_CYLINDER_PF_CACHE_HINT__ 1
#endif 

namespace gms {


          namespace radiolocation {




               namespace {
                   const  static __m512 Ir  = _mm512_setzero_ps();
                   const  static __m512 Ii  = _mm512_set1_ps(1.0f);
                   const  static __m512 nIr = _mm512_set1_ps(-0.0f);
                   const  static __m512 nIi = _mm512_set1_ps(-1.0f);
                   const  static __m512 PI  = _mm512_set1_ps(3.14159265358979323846264338328f);

               }


                   /* 
                         Low frequency scattering widths (k0a << 1).
                         Backscatter scattering width for E-field 
                         cylinder-parallel,formula 4.1-19
                    */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f419_zmm16r4(const __m512 a,
                                           const __m512 k0a) {

                          const register __m512 num = _mm512_mul_ps(a, 
                                                           _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 pi4 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const register __m512 c0  = _mm512_set1_ps(0.8905f);
                          const register __m512 arg = _mm512_mul_ps(k0a,c0);
                          __m512 ln,ln2,rcs,den;
                          ln = logkf(arg);
                          ln2= _mm512_mul_ps(ln,ln);
                          den= _mm512_fmadd_ps(k0a,ln2,pi4);
                          rcs= _mm512_div_ps(num,den);
                          return (rcs);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f419_zmm16r4_unroll16x(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 16;
                         register __m512 a,k0a,rcs;
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
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+0] = rcs;  
                              a   = pa[j+1];
                              k0a = pk0a[j+1];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+1] = rcs;
                              a   = pa[j+2];
                              k0a = pk0a[j+2];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+2] = rcs; 
                              a   = pa[j+3];
                              k0a = pk0a[j+3];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+3] = rcs;   
                              a   = pa[j+4];
                              k0a = pk0a[j+4];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+4] = rcs;
                              a   = pa[j+5];
                              k0a = pk0a[j+5];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+5] = rcs; 
                              a   = pa[j+6];
                              k0a = pk0a[j+6];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+6] = rcs;   
                              a   = pa[j+7];
                              k0a = pk0a[j+7];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+7] = rcs;
                              a   = pa[j+8];
                              k0a = pk0a[j+8];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+8] = rcs;  
                              a   = pa[j+9];
                              k0a = pk0a[j+9];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+9] = rcs;
                              a   = pa[j+10];
                              k0a = pk0a[j+10];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+10] = rcs; 
                              a   = pa[j+11];
                              k0a = pk0a[j+11];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+11] = rcs;
                              a   = pa[j+12];
                              k0a = pk0a[j+12];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+12] = rcs;
                              a   = pa[j+13];
                              k0a = pk0a[j+13];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+13] = rcs;
                              a   = pa[j+14];
                              k0a = pk0a[j+14];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+14] = rcs;
                              a   = pa[j+15];
                              k0a = pk0a[j+15];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+15] = rcs;
                         }            
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f419_zmm16r4_unroll10x(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 10;
                         register __m512 a,k0a,rcs;
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
                         for(j = m1; j != n; j += 10) {[
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
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+0] = rcs;  
                              a   = pa[j+1];
                              k0a = pk0a[j+1];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+1] = rcs;
                              a   = pa[j+2];
                              k0a = pk0a[j+2];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+2] = rcs; 
                              a   = pa[j+3];
                              k0a = pk0a[j+3];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+3] = rcs;   
                              a   = pa[j+4];
                              k0a = pk0a[j+4];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+4] = rcs;
                              a   = pa[j+5];
                              k0a = pk0a[j+5];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+5] = rcs; 
                              a   = pa[j+6];
                              k0a = pk0a[j+6];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+6] = rcs;   
                              a   = pa[j+7];
                              k0a = pk0a[j+7];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+7] = rcs;
                              a   = pa[j+8];
                              k0a = pk0a[j+8];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+8] = rcs;  
                              a   = pa[j+9];
                              k0a = pk0a[j+9];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+9] = rcs;
                           
                         }            
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f419_zmm16r4_unroll6x(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 6;
                         register __m512 a,k0a,rcs;
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
                         for(j = m1; j != n; j += 6) {[
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
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+0] = rcs;  
                              a   = pa[j+1];
                              k0a = pk0a[j+1];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+1] = rcs;
                              a   = pa[j+2];
                              k0a = pk0a[j+2];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+2] = rcs; 
                              a   = pa[j+3];
                              k0a = pk0a[j+3];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+3] = rcs;   
                              a   = pa[j+4];
                              k0a = pk0a[j+4];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+4] = rcs;
                              a   = pa[j+5];
                              k0a = pk0a[j+5];
                              rcs = rcs_f419_zmm16r4(a,k0a);
                              prcs[j+5] = rcs; 
                                                       
                         }            
              }

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f419_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 a   = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const register __m512 num = _mm512_mul_ps(a, 
                                                           _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 pi4 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const register __m512 c0  = _mm512_set1_ps(0.8905f);
                          const register __m512 arg = _mm512_mul_ps(k0a,c0);
                          __m512 ln,ln2,rcs,den;
                          ln = logkf(arg);
                          ln2= _mm512_mul_ps(ln,ln);
                          den= _mm512_fmadd_ps(k0a,ln2,pi4);
                          rcs= _mm512_div_ps(num,den);
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f419_zmm16r4_u(const float * __restrict  pa,
                                             const float * __restrict  pk0a) {

                          const register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 num = _mm512_mul_ps(a, 
                                                           _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 pi4 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const register __m512 c0  = _mm512_set1_ps(0.8905f);
                          const register __m512 arg = _mm512_mul_ps(k0a,c0);
                          __m512 ln,ln2,rcs,den;
                          ln = logkf(arg);
                          ln2= _mm512_mul_ps(ln,ln);
                          den= _mm512_fmadd_ps(k0a,ln2,pi4);
                          rcs= _mm512_div_ps(num,den);
                          return (rcs);
              }


                /* 
                         Low frequency scattering widths (k0a << 1).
                         Backscatter scattering width for H-field 
                         cylinder-parallel,formula 4.1-20
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4120_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          const register __m512 pi2a = _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 c0   = _mm512_set1_ps(2.25f);
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,
                                                                 _mm512_mul_ps(k0a,k0a));
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(c0,k0a3));
                          return (rcs);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f4120_zmm16r4_unroll16x(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST) {
                                                   
                         if(__builtin_expect(n<=0,0)) {return;}
                         if(__builtin_expect(PF_DIST<=0)) PF_DIST = 16;
                         register __m512 a,k0a,rcs;
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
                   __m512 rcs_f4120_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 a    = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          const register __m512 pi2a = _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 c0   = _mm512_set1_ps(2.25f);
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,
                                                                 _mm512_mul_ps(k0a,k0a));
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(c0,k0a3));
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4120_zmm16r4_u(const float * __restrict  pa,
                                            const float * __restrict  pk0a) {

                          const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 pi2a = _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 c0   = _mm512_set1_ps(2.25f);
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,
                                                                 _mm512_mul_ps(k0a,k0a));
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(c0,k0a3));
                          return (rcs);
              }


                /*
                        Bistatic scattering widths, E-field cylinder axis-parallel
                        Formula 4.1-21
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4121_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          return (rcs_f4120_zmm16r4(a,k0a));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4121_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          return (rcs_f4120_zmm16r4_a(pa,pk0a));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4121_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          return (rcs_f4120_zmm16r4_u(pa,pk0a));
              }



                 /*
                        Bistatic scattering widths, H-field cylinder axis-parallel
                        Formula 4.1-22
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4122_zmm16r4(const __m512 phi,
                                            const __m512 a,
                                            const __m512 k0a) {

                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_ps(hlf,cosph);
                          sqr   = _mm512_mul_ps(frac,frac);
                          rcs   = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4122_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 phi = _mm512_load_ps(&pphi[0]);
                          const register __m512 a   = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_ps(hlf,cosph);
                          sqr   = _mm512_mul_ps(frac,frac);
                          rcs   = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4122_zmm16r4_u(const float * __restrict  pphi,
                                              const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          const register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                          const register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_ps(hlf,cosph);
                          sqr   = _mm512_mul_ps(frac,frac);
                          rcs   = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                 }


                   /*
                       Forward scattering widths, E-field.
                       Formula 4.1-23
                   */
 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4123_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          return (rcs_f4120_zmm16r4(a,k0a));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4123_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          return (rcs_f4120_zmm16r4_a(pa,pk0a));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4123_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          return (rcs_f4120_zmm16r4_u(pa,pk0a));
              }


                  /*
                       Forward scattering widths, H-field.
                       Formula 4.1-24
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4124_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const register __m512 qtr = _mm512_set1_ps(0.25f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,qtr));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4124_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 a   = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const register __m512 qtr = _mm512_set1_ps(0.25f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,qtr));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4124_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          const register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const register __m512 qtr = _mm512_set1_ps(0.25f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,qtr));
                          return (rcs);
                  }


                    /*
                          Surface currents (k0a << 1), for long cylinder (wire).
                          E-field cylinder axis parallel.
                          Formula 4.1-25
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm16r4(const __m512 eps0,
                                         const __m512 mu0,
                                         const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 k0a,
                                         __m512 * __restrict Kzr,
                                         __m512 * __restrict Kzi) {

                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0a,ln);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&Kzr,&Kzi);
                      
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) peps0,
                                           const  float * __restrict __ATTR_ALIGN__(64) pmu0,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const   float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) Kzr,
                                           float * __restrict __ATTR_ALIGN__(64) Kzi) {

                        const register __m512 eps0 = _mm512_load_ps(&peps0[0]);
                        const register __m512 mu0  = _mm512_load_ps(&pmu0[0]);
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0a,ln);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_ps(&Kzr[0], resr);
                        _mm512_store_ps(&Kzi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm16r4_u(const  float * __restrict  peps0,
                                           const  float * __restrict  pmu0,
                                           const   float * __restrict  pEr,
                                           const   float * __restrict  pEi,
                                           const   float * __restrict  pk0a,
                                           float * __restrict  Kzr,
                                           float * __restrict  Kzi) {

                        const register __m512 eps0 = _mm512_loadu_ps(&peps0[0]);
                        const register __m512 mu0  = _mm512_loadu_ps(&pmu0[0]);
                        const register __m512 Er   = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_loadu_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0a,ln);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_storeu_ps(&Kzr[0], resr);
                        _mm512_storeu_ps(&Kzi[0], resi);
                 }


                  /*
                          Surface currents (k0a << 1), for long cylinder (wire).
                          H-field cylinder axis parallel.
                          Formula 4.1-26
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm16r4(const __m512 Hr,
                                          const __m512 Hi,
                                          __m512 * __restrict Kphr,
                                          __m512 * __restrict Kphi) {

                        *Kphr = _mm512_mul_ps(nIi,Hr);
                        *Kphi = _mm512_mul_ps(nIi,Hi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) Hr,
                                            const float * __restrict __ATTR_ALIGN__(64) Hi,
                                           float * __restrict __ATTR_ALIGN__(64) Kphr,
                                          float * __restrict __ATTR_ALIGN__(64) Kphi) {

                        _mm512_store_ps(&Kphr[0] ,_mm512_mul_ps(nIi,_mm512_load_ps(&Hr[0]));
                        _mm512_store_ps(&Kphi[0] ,_mm512_mul_ps(nIi,_mm512_load_ps(&Hi[0]));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm16r4_u(const float * __restrict  Hr,
                                            const float * __restrict  Hi,
                                           float * __restrict  Kphr,
                                          float * __restrict Kphi) {

                        _mm512_storeu_ps(&Kphr[0] ,_mm512_mul_ps(nIi,_mm512_loadu_ps(&Hr[0]));
                        _mm512_storeu_ps(&Kphi[0] ,_mm512_mul_ps(nIi,_mm512_loadu_ps(&Hi[0]));
                 }


                   /*
                        The toal current along the wire.
                        Formula 4.1-27 

                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm16r4(const __m512 eps0,
                                         const __m512 mu0,
                                         const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 k0a,
                                         const __m512 k0,
                                         __m512 * __restrict Izr,
                                         __m512 * __restrict Izi) {

                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 _2pi= _mm512_set1_ps(6.283185307179586476925286766559f);
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t2r = _mm512_mul_ps(_2pi,Er);
                        t1r = _mm512_mul_ps(k0,ln);
                        t2i = _mm512_mul_ps(_2pi,Ei);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&Izr,&Izi);
                      
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) peps0,
                                           const  float * __restrict __ATTR_ALIGN__(64) pmu0,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const   float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const   float * __restrict __ATTR_ALIGN__(64) pk0,
                                           float * __restrict __ATTR_ALIGN__(64) Izr,
                                           float * __restrict __ATTR_ALIGN__(64) Izi) {

                        const register __m512 eps0 = _mm512_load_ps(&peps0[0]);
                        const register __m512 _2pi= _mm512_set1_ps(6.283185307179586476925286766559f);
                        const register __m512 mu0  = _mm512_load_ps(&pmu0[0]);
                        const register __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 ln  = logkf(lna);
                        const register __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0,ln);
                        t2r = _mm512_mul_ps(_2pi,Er);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        t2i = _mm512_mul_ps(_2pi,Ei);
                        cdiv_zmm16r4(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_ps(&Izr[0], resr);
                        _mm512_store_ps(&Izi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm16r4_u(const  float * __restrict  peps0,
                                           const  float * __restrict  pmu0,
                                           const   float * __restrict  pEr,
                                           const   float * __restrict  pEi,
                                           const   float * __restrict  pk0a,
                                           const   float * __restrict  pk0,
                                           float * __restrict  Izr,
                                           float * __restrict  Izi) {

                        const register __m512 eps0 = _mm512_load_ps(&peps0[0]);
                        const register __m512 _2pi= _mm512_set1_ps(6.283185307179586476925286766559f);
                        const register __m512 mu0  = _mm512_load_ps(&pmu0[0]);
                        const register __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 ln  = logkf(lna);
                        const register __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0,ln);
                        t2r = _mm512_mul_ps(_2pi,Er);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        t2i = _mm512_mul_ps(_2pi,Ei);
                        cdiv_zmm16r4(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_ps(&Izr[0], resr);
                        _mm512_store_ps(&Izi[0], resi);
                 }


                   /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        (0<<phi<pi/2, k0a > 2)
                        Electric-field.
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm16r4(const __m512 phi2,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 Er,
                                         const __m512 Ei,
                                         __m512 * __restrict EOr,
                                         __m512 * __restrict EOi) {

                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,exr,exi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(5.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_sub_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_ps(c6,t2i);
                        cmul_zmm16r4(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_ps(t4,_mm512_add_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(exr,exi,t2r,t2i,&t0r,&t0i);
                        *EOr = t0r;
                        *EOi = t0i;
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           float * __restrict __ATTR_ALIGN__(64) EOr,
                                           float * __restrict __ATTR_ALIGN__(64) EOi) {

                        const register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,exr,exi;;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(5.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_sub_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_ps(c6,t2i);
                        cmul_zmm16r4(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_ps(t4,_mm512_add_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(exr,exi,t2r,t2i,&t0r,&t0i);
                        _mm512_store_ps(&EOr[0], t0r);
                        _mm512_store_ps(&EOi[0], t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm16r4_u(const float * __restrict  pphi2,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           float * __restrict  EOr,
                                           float * __restrict  EOi) {

                        const register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_loadu_ps(&pEi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,exr,exi;;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(5.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_sub_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_ps(c6,t2i);
                        cmul_zmm16r4(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_ps(t4,_mm512_add_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(exr,exi,t2r,t2i,&t0r,&t0i);
                        _mm512_storeu_ps(&EOr[0], t0r);
                        _mm512_storeu_ps(&EOi[0], t0i);
                 }


                     /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        (0<<phi<pi/2, k0a > 2)
                        Magnetic-field.
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm16r4(const __m512 phi2,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 Hr,
                                         const __m512 Hi,
                                         __m512 * __restrict HOr,
                                         __m512 * __restrict HOi) {

                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(7.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_add_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_ps(c6,t2i);
                        cmul_zmm16r4(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_ps(t4,_mm512_sub_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        *HOr = t0r;
                        *HOi = t0i;
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pHr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHi,
                                           float * __restrict __ATTR_ALIGN__(64) HOr,
                                           float * __restrict __ATTR_ALIGN__(64) HOi) {

                        const register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_load_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_load_ps(&pHi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(7.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_add_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_ps(c6,t2i);
                        cmul_zmm16r4(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_ps(t4,_mm512_sub_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        _mm512_store_ps(&HOr[0], t0r);
                        _mm512_store_ps(&HOi[0], t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm16r4_u(const float * __restrict  pphi2,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pHr,
                                           const float * __restrict  pHi,
                                           float * __restrict  HOr,
                                           float * __restrict  HOi) {

                        const register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_loadu_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_loadu_ps(&pHi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(7.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_add_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_ps(c6,t2i);
                        cmul_zmm16r4(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_ps(t4,_mm512_sub_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        _mm512_storeu_ps(&HOr[0], t0r);
                        _mm512_storeu_ps(&HOi[0], t0i);
                 }


                 
                 /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        (0<<phi<pi/2, k0a > 2)
                        Electric-field.
                        Formula 4.1-30
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         __m512 * __restrict ECr,
                                         __m512 * __restrict ECi) {

                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(0.910721f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.9358135f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(1.607129f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.057397f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.0994145f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        *ECr = _mm512_add_ps(t0r,tmp3r);
                        *ECi = _mm512_add_ps(t0i,tmp3i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi,
                                           float * __restrict __ATTR_ALIGN__(64) ECr,
                                           float * __restrict __ATTR_ALIGN__(64) ECi) {

                        const register __m512 phi  = _mm512_load_ps(&pphi[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(0.910721f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.9358135f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(1.607129f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.057397f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.0994145f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_store_ps(&ECr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_store_ps(&ECi[0] ,_mm512_add_ps(t0i,tmp3i));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm16r4_a(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pphi,
                                           float * __restrict  ECr,
                                           float * __restrict ECi) {

                        const register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_loadu_ps(&pEi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(0.910721f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.9358135f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(1.607129f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.057397f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.0994145f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_storeu_ps(&ECr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_storeu_ps(&ECi[0] ,_mm512_add_ps(t0i,tmp3i));
                 }


                    /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        valid only for (0<<phi<pi/2, k0a > 2)
                        Magnetic-field.
                        Formula 4.1-32
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm16r4(const __m512 Hr,
                                         const __m512 Hi,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         __m512 * __restrict HCr,
                                         __m512 * __restrict HCi) {

                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(1.531915f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.404308f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(0.70028f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.072732f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.1259755f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        *HCr = _mm512_add_ps(t0r,tmp3r);
                        *HCi = _mm512_add_ps(t0i,tmp3i);
                 }


                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64)  pHr,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pHi,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pa,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           float * __restrict __ATTR_ALIGN__(64)  HCr,
                                           float * __restrict __ATTR_ALIGN__(64)  HCi) {

                        const register __m512 phi  = _mm512_load_ps(&pphi[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_load_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_load_ps(&pHi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(1.531915f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.404308f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(0.70028f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.072732f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.1259755f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_store_ps(&HCr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_store_ps(&HCi[0] ,_mm512_add_ps(t0i,tmp3i));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm16r4_u(const  float * __restrict  pHr,
                                           const  float * __restrict  pHi,
                                           const  float * __restrict  pa,
                                           const  float * __restrict  pr,
                                           const  float * __restrict  pk0,
                                           const  float * __restrict  pk0a,
                                           const  float * __restrict  pphi,
                                           float * __restrict   HCr,
                                           float * __restrict   HCi) {

                        const register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_loadu_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_loadu_ps(&pHi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(1.531915f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.404308f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(0.70028f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.072732f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.1259755f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_storeu_ps(&HCr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_storeu_ps(&HCi[0] ,_mm512_add_ps(t0i,tmp3i));
                 }


                   /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Optical wave component e-field, formula 4.1-33
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0a,
                                         __m512 * __restrict EOr,
                                         __m512 * __restrict EOi) {

                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _5 = _mm512_set1_ps(5.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(127.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,*EOr,*EOi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) EOr,
                                           float * __restrict __ATTR_ALIGN__(64) EOi) {

                        register __m512 Er = _mm512_load_ps(&pEr[0]);
                        register __m512 Ei = _mm512_load_ps(&pEi[0]);
                        register __m512 a  = _mm512_load_ps(&pa[0]);
                        register __m512 r  = _mm512_load_ps(&pr[0]);
                        register __m512 k0a= _mm512_load_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        register __m512 resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _5 = _mm512_set1_ps(5.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(127.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&EOr[0],resr);
                        _mm512_store_ps(&EOi[0],resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm16r4_u(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0a,
                                           float * __restrict  EOr,
                                           float * __restrict  EOi) {

                        register __m512 Er = _mm512_loadu_ps(&pEr[0]);
                        register __m512 Ei = _mm512_loadu_ps(&pEi[0]);
                        register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0a= _mm512_loadu_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        register __m512 resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _5 = _mm512_set1_ps(5.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(127.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&EOr[0],resr);
                        _mm512_storeu_ps(&EOi[0],resi);
                 }


                     /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Optical wave component h-field, formula 4.1-35
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm16r4(const __m512 Hr,
                                         const __m512 Hi,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0a,
                                         __m512 * __restrict HOr,
                                         __m512 * __restrict HOi) {

                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _11 = _mm512_set1_ps(11.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(353.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,*HOr,*HOi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pHr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) HOr,
                                           float * __restrict __ATTR_ALIGN__(64) HOi) {

                        register __m512 Hr = _mm512_load_ps(&pHr[0]);
                        register __m512 Hi = _mm512_load_ps(&pHi[0]);
                        register __m512 a  = _mm512_load_ps(&pa[0]);
                        register __m512 r  = _mm512_load_ps(&pr[0]);
                        register __m512 k0a= _mm512_load_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i,resr,resi;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _11 = _mm512_set1_ps(11.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(353.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&HOr[0], resr);
                        _mm512_store_ps(&HOi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm16r4_u(const float * __restrict  pHr,
                                           const float * __restrict  pHi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0a,
                                           float * __restrict  HOr,
                                           float * __restrict  HOi) {

                        register __m512 Hr = _mm512_loadu_ps(&pHr[0]);
                        register __m512 Hi = _mm512_loadu_ps(&pHi[0]);
                        register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0a= _mm512_loadu_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i,resr,resi;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _11 = _mm512_set1_ps(11.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(353.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&HOr[0], resr);
                        _mm512_storeu_ps(&HOi[0], resi);
                 }


                   /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Creeping wave component e-field, formula 4.1-34
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         __m512 * __restrict ECr,
                                         __m512 * __restrict ECi) {

                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(2.939945f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.180318f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(1.821442f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-5.048945f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.312320f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Er,t2);
                        fraci = _mm512_mul_ps(Ei,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        *ECr = _mm512_mul_ps(t0r,rex);
                        *ECi = _mm512_mul_ps(t0i,rex);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                         const float * __restrict __ATTR_ALIGN__(64) pEi,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pr,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         float * __restrict __ATTR_ALIGN__(64) ECr,
                                         float * __restrict __ATTR_ALIGN__(64) ECi) {

                        const register __m512 Er = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei = _mm512_load_ps(&pEi[0]);
                        const register __m512 a  = _mm512_load_ps(&pa[0]);
                        const register __m512 r  = _mm512_load_ps(&pr[0]);
                        const register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(2.939945f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.180318f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(1.821442f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-5.048945f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.312320f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Er,t2);
                        fraci = _mm512_mul_ps(Ei,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_store_ps(&ECr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_store_ps(&ECi[0] ,_mm512_mul_ps(t0i,rex));
                }


                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm16r4_u(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           float * __restrict  ECr,
                                           float * __restrict  ECi) {

                        const register __m512 Er = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei = _mm512_loadu_ps(&pEi[0]);
                        const register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(2.939945f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.180318f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(1.821442f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-5.048945f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.312320f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Er,t2);
                        fraci = _mm512_mul_ps(Ei,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_storeu_ps(&ECr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_storeu_ps(&ECi[0] ,_mm512_mul_ps(t0i,rex));
                }



                   /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Creeping wave component h-field, formula 4.1-36
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm16r4(const __m512 Hr,
                                         const __m512 Hi,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         __m512 * __restrict HCr,
                                         __m512 * __restrict HCi) {

                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(1.2701695f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.2284945f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(3.063830f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-2.200000f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.3957635f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Hr,t2);
                        fraci = _mm512_mul_ps(Hi,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        *HCr = _mm512_mul_ps(t0r,rex);
                        *HCi = _mm512_mul_ps(t0i,rex);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pHr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           float * __restrict __ATTR_ALIGN__(64) HCr,
                                           float * __restrict __ATTR_ALIGN__(64) HCi) {

                        const register __m512 Hr = _mm512_load_ps(&pHr[0]);
                        const register __m512 Hi = _mm512_load_ps(&pHi[0]);
                        const register __m512 a  = _mm512_load_ps(&pa[0]);
                        const register __m512 r  = _mm512_load_ps(&pr[0]);
                        const register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(1.2701695f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.2284945f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(3.063830f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-2.200000f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.3957635f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Hr,t2);
                        fraci = _mm512_mul_ps(Hi,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_store_ps(&HCr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_store_ps(&HCi[0] ,_mm512_mul_ps(t0i,rex));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm16r4_u(const float * __restrict  pHr,
                                           const float * __restrict  pHi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           float * __restrict  HCr,
                                           float * __restrict  HCi) {

                        const register __m512 Hr = _mm512_loadu_ps(&pHr[0]);
                        const register __m512 Hi = _mm512_loadu_ps(&pHi[0]);
                        const register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(1.2701695f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.2284945f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(3.063830f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-2.200000f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.3957635f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Hr,t2);
                        fraci = _mm512_mul_ps(Hi,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_storeu_ps(&HCr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_storeu_ps(&HCi[0] ,_mm512_mul_ps(t0i,rex));
                }


                  /*
                        Bistatic scattering width in high frequency limit (k0a > 20)
                        for |PI-phi| > k0a^0.3
                        Formula 4.1-37
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4137_zmm16r4(const __m512 a,
                                            const __m512 phi2) {

                          register __m512 rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_ps(PI,_mm512_mul_ps(a,cosp2));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4137_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi2) {

                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                          register __m512 rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_ps(PI,_mm512_mul_ps(a,cosp2));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4137_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pphi2) {

                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                          register __m512 rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_ps(PI,_mm512_mul_ps(a,cosp2));
                          return (rcs);
                 }


                    /*
                          Backscattering Width in High-Frequency Limit (k0a > 20)
                          Formula 4.1-38
                     */

                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4138_zmm16r4(const __m512 a) {

                          return (__m512_mul_ps(a,PI));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4138_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512 a = _mm512_load_ps(&pa[0]);
                          return (__m512_mul_ps(a,PI));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4138_zmm16r4_u(const float * __restrict  pa) {

                          register __m512 a = _mm512_loadu_ps(&pa[0]);
                          return (__m512_mul_ps(a,PI));
                  }


                   /*
                         Forward scattering widths and pattern in high-frequency limit
                         (k0a>20.0)
                         Formula 4.1-40, RCS.
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4140_zmm16r4(const __m512 k0a,
                                            const __m512 alpha) {

                          register __m512 sinc,k0alp,k0as,t0;
                          register __m512 rcs;
                          const __m512 _4       = _mm512_set1_ps(4.0f);
                          k0alp = _mm512_mul_ps(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_ps(t0,k0alp); 
                          k0as  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a)); 
                          rcs   = _mm512_mul_ps(k0as,_mm512_mul_ps(sinc,sinc));
                          return (rcs);        
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4140_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) palpha) {

                          register __m512 k0a   = _mm512_load_ps(&pk0a[0]);
                          register __m512 alpha = _mm512_load_ps(&palpha[0]);
                          register __m512 sinc,k0alp,k0as,t0;
                          register __m512 rcs;
                          const __m512 _4       = _mm512_set1_ps(4.0f);
                          k0alp = _mm512_mul_ps(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_ps(t0,k0alp); 
                          k0as  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a)); 
                          rcs   = _mm512_mul_ps(k0as,_mm512_mul_ps(sinc,sinc));
                          return (rcs);        
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4140_zmm16r4_u(const float * __restrict  pk0a,
                                              const float * __restrict  palpha) {

                          register __m512 k0a   = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 alpha = _mm512_loadu_ps(&palpha[0]);
                          register __m512 sinc,k0alp,k0as,t0;
                          register __m512 rcs;
                          const __m512 _4       = _mm512_set1_ps(4.0f);
                          k0alp = _mm512_mul_ps(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_ps(t0,k0alp); 
                          k0as  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a)); 
                          rcs   = _mm512_mul_ps(k0as,_mm512_mul_ps(sinc,sinc));
                          return (rcs);        
                  }


                     /*
                         Forward scattering widths and pattern in high-frequency limit
                         (k0a>20.0), forward scattered (diffracted) e-field
                         Formula 4.1-39.

                       */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 alp
                                         const __m512 k0a,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                        register __m512 _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512 facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512 cer,cei;
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        _2k0as = _mm512_add_ps(_2,_mm512_mul_ps(k0a,k0a));
                        k0r    = _mm512_mul_ps(k0,r);
                        k0alp  = _mm512_mul_ps(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_ps(k0r,pir);
                        sinc   = _mm512_div_ps(t0,k0alp);
                        cexp_zmm16r4(arr,ari,&cer,&cei);
                        div    = _mm512_div_ps(_2k0as,pi4);
                        t0     = _mm512_sqrt_ps(div);
                        facr   = _mm512_mul_ps(Er,t0);
                        t0r    = _mm512_mul_ps(cer,sinc);
                        faci   = _mm512_mul_ps(Ei,t0);
                        t0i    = _mm512_mul_ps(cei,sinc);
                        cmul_zmm16r4(facr,faci,t0r,t0i,*Esr,*Esi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) palp
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) Esr,
                                           float * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512 Er = _mm512_load_ps(&pEr[0]);
                        register __m512 Ei = _mm512_load_ps(&pEi[0]);
                        register __m512 r  = _mm512_load_ps(&pr[0]);
                        register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 alp= _mm512_load_ps(&palp[0]);
                        register __m512 k0a= _mm512_load_ps(&pk0a[0]);
                        register __m512 _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512 facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512 cer,cei,resr,resi;
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        _2k0as = _mm512_add_ps(_2,_mm512_mul_ps(k0a,k0a));
                        k0r    = _mm512_mul_ps(k0,r);
                        k0alp  = _mm512_mul_ps(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_ps(k0r,pir);
                        sinc   = _mm512_div_ps(t0,k0alp);
                        cexp_zmm16r4(arr,ari,&cer,&cei);
                        div    = _mm512_div_ps(_2k0as,pi4);
                        t0     = _mm512_sqrt_ps(div);
                        facr   = _mm512_mul_ps(Er,t0);
                        t0r    = _mm512_mul_ps(cer,sinc);
                        faci   = _mm512_mul_ps(Ei,t0);
                        t0i    = _mm512_mul_ps(cei,sinc);
                        cmul_zmm16r4(facr,faci,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&Esr[0], resr);
                        _mm512_store_ps(&Esi[0], resi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm16r4_u(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  palp
                                           const float * __restrict  pk0a,
                                           float * __restrict  Esr,
                                           float * __restrict  Esi) {

                        register __m512 Er = _mm512_loadu_ps(&pEr[0]);
                        register __m512 Ei = _mm512_loadu_ps(&pEi[0]);
                        register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 alp= _mm512_loadu_ps(&palp[0]);
                        register __m512 k0a= _mm512_loadu_ps(&pk0a[0]);
                        register __m512 _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512 facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512 cer,cei,resr,resi;
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        _2k0as = _mm512_add_ps(_2,_mm512_mul_ps(k0a,k0a));
                        k0r    = _mm512_mul_ps(k0,r);
                        k0alp  = _mm512_mul_ps(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_ps(k0r,pir);
                        sinc   = _mm512_div_ps(t0,k0alp);
                        cexp_zmm16r4(arr,ari,&cer,&cei);
                        div    = _mm512_div_ps(_2k0as,pi4);
                        t0     = _mm512_sqrt_ps(div);
                        facr   = _mm512_mul_ps(Er,t0);
                        t0r    = _mm512_mul_ps(cer,sinc);
                        faci   = _mm512_mul_ps(Ei,t0);
                        t0i    = _mm512_mul_ps(cei,sinc);
                        cmul_zmm16r4(facr,faci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&Esr[0], resr);
                        _mm512_storeu_ps(&Esi[0], resi);
              }


                  /*
                         Forward scattering widths and pattern in high-frequency limit
                         (k0a>20.0), constant angle (alpha=0)
                         Formula 4.1-41, RCS.
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4141_zmm16r4(const __m512 k0a) {

                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a));
                          return (rcs);
                 }

                

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4141_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0a) {

                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4141_zmm16r4_u(const float * __restrict   pk0a) {

                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a));
                          return (rcs);
                 }


                   /*
                        Approximations for the low frequency region (k0a<<1,k1a<<1)
                        Scattered far-zone e-field, formula 4.1-45
                    */

                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm16r4(const __m512 EIr,
                                         const __m512 EIi,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         const __m512 eps0,
                                         const __m512 eps1,
                                         const __m512 mu0,
                                         const __m512 mu1,
                                         __m512 * __restrict ESr,
                                         __m512 * __restrict ESi) {

                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(EIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                        fraci           = _mm512_mul_ps(EIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(mu1,mu0),
                                                        _mm512_add_ps(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,*ESr,*ESi);
               }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEIr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps0,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps1,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu1,
                                           float * __restrict __ATTR_ALIGN__(64)  ESr,
                                           float * __restrict __ATTR_ALIGN__(64)  ESi) {

                        register __m512 EIr = _mm512_load_ps(&pEIr[0]);
                        register __m512 EIi = _mm512_load_ps(&pEIi[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 eps0= _mm512_load_ps(&peps0[0]);
                        register __m512 eps1= _mm512_load_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(EIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                        fraci           = _mm512_mul_ps(EIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(mu1,mu0),
                                                        _mm512_add_ps(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&ESr[0], resr);
                        _mm512_store_ps(&ESi[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm16r4_u(const float * __restrict  pEIr,
                                           const float * __restrict  pEIi,
                                           const float * __restrict  pr,
                                           const float * __restrict   pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict   pphi,
                                           const float * __restrict   peps0,
                                           const float * __restrict  peps1,
                                           const float * __restrict   pmu0,
                                           const float * __restrict   pmu1,
                                           float * __restrict   ESr,
                                           float * __restrict   ESi) {

                        register __m512 EIr = _mm512_loadu_ps(&pEIr[0]);
                        register __m512 EIi = _mm512_loadu_ps(&pEIi[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                        register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(EIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                        fraci           = _mm512_mul_ps(EIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(mu1,mu0),
                                                        _mm512_add_ps(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&ESr[0], resr);
                        _mm512_storeu_ps(&ESi[0], resi);
               }


                  /*
                        Approximations for the low frequency region (k0a<<1,k1a<<1)
                        Scattered far-zone h-field, formula 4.1-46
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm16r4(const __m512 HIr,
                                         const __m512 HIi,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         const __m512 eps0,
                                         const __m512 eps1,
                                         const __m512 mu0,
                                         const __m512 mu1,
                                         __m512 * __restrict HSr,
                                         __m512 * __restrict HSi) {

                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(HIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                        fraci           = _mm512_mul_ps(HIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(eps1,eps0),
                                                        _mm512_add_ps(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,*HSr,*HSi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pHIr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps0,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps1,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu1,
                                           float * __restrict __ATTR_ALIGN__(64)  HSr,
                                           float * __restrict __ATTR_ALIGN__(64)  HSi) {

                        register __m512 HIr = _mm512_load_ps(&pHIr[0]);
                        register __m512 HIi = _mm512_load_ps(&pHIi[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 eps0= _mm512_load_ps(&peps0[0]);
                        register __m512 eps1= _mm512_load_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(HIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                        fraci           = _mm512_mul_ps(HIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(eps1,eps0),
                                                        _mm512_add_ps(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&HSr[0], resr);
                        _mm512_store_ps(&HSi[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm16r4_u(const float * __restrict   pHIr,
                                           const float * __restrict   pHIi,
                                           const float * __restrict   pr,
                                           const float * __restrict   pk0,
                                           const float * __restrict   pk0a,
                                           const float * __restrict   pphi,
                                           const float * __restrict   peps0,
                                           const float * __restrict   peps1,
                                           const float * __restrict   pmu0,
                                           const float * __restrict   pmu1,
                                           float * __restrict  HSr,
                                           float * __restrict   HSi) {

                        register __m512 HIr = _mm512_loadu_ps(&pHIr[0]);
                        register __m512 HIi = _mm512_loadu_ps(&pHIi[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                        register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(HIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                        fraci           = _mm512_mul_ps(HIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(eps1,eps0),
                                                        _mm512_add_ps(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&HSr[0], resr);
                        _mm512_storeu_ps(&HSi[0], resi);
               }


                 /*
                      Bistatic scattering width (k0a<<1, k1a<<1) at the angle 'phi'
                      Formula 4.1-47

                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4147_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 phi,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }


                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4147_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(64) pphi,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps0,
                                            const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                            const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          register __m512 phi = _mm512_load_ps(&pphi[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 pmu1= _mm512_load_ps(&pmu1[0]);
                          register __m512 pmu0= _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4147_zmm16r4_u(const float * __restrict  pa,
                                            const float * __restrict  pk0a,
                                            const float * __restrict  pphi,
                                            const float * __restrict  peps1,
                                            const float * __restrict  peps0,
                                            const float * __restrict  pmu1,
                                            const float * __restrict  pmu0) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                          register __m512 pmu1= _mm512_loadu_ps(&pmu1[0]);
                          register __m512 pmu0= _mm512_loadu_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }   


                    /*
                      Bistatic scattering width (k0a<<1, k1a<<1) at the angle 'phi'
                      Formula 4.1-48

                   */   

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4148_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 phi,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }  


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4148_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(64) pphi,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps0,
                                            const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                            const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          register __m512 phi = _mm512_load_ps(&pphi[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 pmu1= _mm512_load_ps(&pmu1[0]);
                          register __m512 pmu0= _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4148_zmm16r4_u(const float * __restrict  pa,
                                            const float * __restrict pk0a,
                                            const float * __restrict  pphi,
                                            const float * __restrict  peps1,
                                            const float * __restrict  peps0,
                                            const float * __restrict  pmu1,
                                            const float * __restrict  pmu0) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                          register __m512 pmu1= _mm512_loadu_ps(&pmu1[0]);
                          register __m512 pmu0= _mm512_loadu_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }  


                   /*
                         Backscattering width (k0a<<1,k1a<<1), when phi = 0
                         Formula 4.1-49
                    */
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4149_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4149_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) peps1,
                                              const float * __restrict __ATTR_ALIGN__(64) peps0,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4149_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  peps1,
                                              const float * __restrict  peps0,
                                              const float * __restrict  pmu1,
                                              const float * __restrict  pmu0) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                     /*
                         Backscattering width (k0a<<1,k1a<<1), when phi = 0
                         Formula 4.1-50
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4150_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4150_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) peps1,
                                              const float * __restrict __ATTR_ALIGN__(64) peps0,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4150_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  peps1,
                                              const float * __restrict  peps0,
                                              const float * __restrict  pmu1,
                                              const float * __restrict  pmu0) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                    /*
                         Forward scattering width (k0a<<1, k1a<<1), phi = pi
                         Formula 4.1-51
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4151_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_add_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4151_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) peps1,
                                              const float * __restrict __ATTR_ALIGN__(64) peps0,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_add_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4151_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  peps1,
                                              const float * __restrict  peps0,
                                              const float * __restrict  pmu1,
                                              const float * __restrict  pmu0) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_add_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                      /*
                         Forward scattering width (k0a<<1, k1a<<1), phi = pi
                         Formula 4.1-52
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4152_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_add_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4152_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) peps1,
                                              const float * __restrict __ATTR_ALIGN__(64) peps0,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_add_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4152_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  peps1,
                                              const float * __restrict  peps0,
                                              const float * __restrict  pmu1,
                                              const float * __restrict  pmu0) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a3[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                          register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                          register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512 rcs;
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                          t1               = _mm512_sub_ps(eps1,eps0);
                          t2               = _mm512_add_ps(eps1,eps0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_add_ps(epst,mut);
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);

                   }


                     /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-72
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4172_zmm16r4(const __m512 mur,
                                          const __m512 mui,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 psi,
                                          __m512 * __restrict Tinr,
                                          __m512 * __restrict Tini) {

                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        register __m512 sin2p,cosp,divr,divi,t1;
                        register __m512 sqr1,sqi1,sqr2,sqi2,t0;
                        register __m512 mulr,muli,t0r,t0i,t1r,t1i;
                        register __m512 t2r,t2i,t3r,t3i;
                        cosp = xcosf(psi);
                        t0   = xsinf(psi);
                        sin2p= _mm512_mul_ps(t0,t0);
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        t1   = _mm512_sub_ps(_1,sin2p);
                        cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                        csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                        t0r = _mm512_div_ps(t1,mulr);
                        t0i = _mm512_div_ps(t1,muli);
                        csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                        t2r = _mm512_add_ps(sqr1,sqr1);
                        t2i = _mm512_add_ps(sqi1,sqi1);
                        cmul_zmm16r4(t2r,t2i,sqr2,sqi2,&t1r,&t1i);//numerator
                        cmul_zmm16r4(sqr1,sqi1,sqr2,sqi2,&t3r,&t3i); // denum
                        t3r = _mm512_add_ps(cosp,t3r);
                        t3i = _mm512_add_ps(cosp,t3i);
                        cdiv_zmm16r4(t1r,t1i,t3r,t3i,*Tinr,*Tini);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4172_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pmur,
                                            const float * __restrict __ATTR_ALIGN__(64) pmui,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                            float * __restrict __ATTR_ALIGN__(64) Tinr,
                                            float * __restrict __ATTR_ALIGN__(64) Tini) {

                        register __m512 mur  = _mm512_load_ps(&pmur[0]);
                        register __m512 mui  = _mm512_load_ps(&pmui[0]);
                        register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        register __m512 sin2p,cosp,divr,divi,t1;
                        register __m512 sqr1,sqi1,sqr2,sqi2,t0;
                        register __m512 mulr,muli,t0r,t0i,t1r,t1i;
                        register __m512 t2r,t2i,t3r,t3i,resr,resi;
                        cosp = xcosf(psi);
                        t0   = xsinf(psi);
                        sin2p= _mm512_mul_ps(t0,t0);
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        t1   = _mm512_sub_ps(_1,sin2p);
                        cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                        csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                        t0r = _mm512_div_ps(t1,mulr);
                        t0i = _mm512_div_ps(t1,muli);
                        csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                        t2r = _mm512_add_ps(sqr1,sqr1);
                        t2i = _mm512_add_ps(sqi1,sqi1);
                        cmul_zmm16r4(t2r,t2i,sqr2,sqi2,&t1r,&t1i);//numerator
                        cmul_zmm16r4(sqr1,sqi1,sqr2,sqi2,&t3r,&t3i); // denum
                        t3r = _mm512_add_ps(cosp,t3r);
                        t3i = _mm512_add_ps(cosp,t3i);
                        cdiv_zmm16r4(t1r,t1i,t3r,t3i,&resr,&resi);
                        _mm512_store_ps(&Tinr[0], resr);
                        _mm512_store_ps(&Tini[0], resi);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4172_zmm16r4_u(const float * __restrict  pmur,
                                            const float * __restrict  pmui,
                                            const float * __restrict  pepsr,
                                            const float * __restrict  pepsi,
                                            const float * __restrict  ppsi,
                                            float * __restrict  Tinr,
                                            float * __restrict  Tini) {

                        register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                        register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        register __m512 sin2p,cosp,divr,divi,t1;
                        register __m512 sqr1,sqi1,sqr2,sqi2,t0;
                        register __m512 mulr,muli,t0r,t0i,t1r,t1i;
                        register __m512 t2r,t2i,t3r,t3i,resr,resi;
                        cosp = xcosf(psi);
                        t0   = xsinf(psi);
                        sin2p= _mm512_mul_ps(t0,t0);
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        t1   = _mm512_sub_ps(_1,sin2p);
                        cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                        csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                        t0r = _mm512_div_ps(t1,mulr);
                        t0i = _mm512_div_ps(t1,muli);
                        csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                        t2r = _mm512_add_ps(sqr1,sqr1);
                        t2i = _mm512_add_ps(sqi1,sqi1);
                        cmul_zmm16r4(t2r,t2i,sqr2,sqi2,&t1r,&t1i);//numerator
                        cmul_zmm16r4(sqr1,sqi1,sqr2,sqi2,&t3r,&t3i); // denum
                        t3r = _mm512_add_ps(cosp,t3r);
                        t3i = _mm512_add_ps(cosp,t3i);
                        cdiv_zmm16r4(t1r,t1i,t3r,t3i,&resr,&resi);
                        _mm512_storeu_ps(&Tinr[0], resr);
                        _mm512_storeu_ps(&Tini[0], resi);
                  }


                     /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-73
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4173_zmm16r4(const __m512 mur,
                                          const __m512 mui,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 psi,
                                          __m512 * __restrict Tinr,
                                          __m512 * __restrict Tini) {

                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         register __m512 cosp,_2cosp,divr,divi;
                         register __m512 sqr1,sqi1,sqr2,sqi2;
                         register __m512 sinp,sin2p,mulr,muli;
                         register __m512 t0r,t0i,_1msp;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         _2cosp = _mm512_add_ps(cosp,cosp);
                         sin2p  = _mm512_mul_ps(sinp,sinp);
                         _1msp  = _mm512_sub_ps(_1,sin2p);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&divi);
                         t0r = _mm512_div_ps(_1msp,mulr);
                         t0i = _mm512_div_ps(_1msp,muli);
                         csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                         csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                         *Tinr = _mm512_fmadd_ps(sqr1,sqr2,cosp);
                         *Tini = _mm512_fmadd_ps(sqi1,sqi2,cosp);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4173_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pmur,
                                            const float * __restrict __ATTR_ALIGN__(64) pmui,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                            float * __restrict __ATTR_ALIGN__(64) Tinr,
                                            float * __restrict __ATTR_ALIGN__(64) Tini) {

                         register __m512 mur  = _mm512_load_ps(&pmur[0]);
                         register __m512 mui  = _mm512_load_ps(&pmui[0]);
                         register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         register __m512 cosp,_2cosp,divr,divi;
                         register __m512 sqr1,sqi1,sqr2,sqi2;
                         register __m512 sinp,sin2p,mulr,muli;
                         register __m512 t0r,t0i,_1msp;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         _2cosp = _mm512_add_ps(cosp,cosp);
                         sin2p  = _mm512_mul_ps(sinp,sinp);
                         _1msp  = _mm512_sub_ps(_1,sin2p);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&divi);
                         t0r = _mm512_div_ps(_1msp,mulr);
                         t0i = _mm512_div_ps(_1msp,muli);
                         csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                         csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                         _mm512_store_ps(&Tinr[0] ,_mm512_fmadd_ps(sqr1,sqr2,cosp));
                         _mm512_store_ps(&Tini[0] ,_mm512_fmadd_ps(sqi1,sqi2,cosp));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4173_zmm16r4_u(const float * __restrict  pmur,
                                            const float * __restrict  pmui,
                                            const float * __restrict  pepsr,
                                            const float * __restrict  pepsi,
                                            const float * __restrict  ppsi,
                                            float * __restrict  Tinr,
                                            float * __restrict  Tini) {

                         register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                         register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         register __m512 cosp,_2cosp,divr,divi;
                         register __m512 sqr1,sqi1,sqr2,sqi2;
                         register __m512 sinp,sin2p,mulr,muli;
                         register __m512 t0r,t0i,_1msp;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         _2cosp = _mm512_add_ps(cosp,cosp);
                         sin2p  = _mm512_mul_ps(sinp,sinp);
                         _1msp  = _mm512_sub_ps(_1,sin2p);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&divi);
                         t0r = _mm512_div_ps(_1msp,mulr);
                         t0i = _mm512_div_ps(_1msp,muli);
                         csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                         csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                         _mm512_storeu_ps(&Tinr[0] ,_mm512_fmadd_ps(sqr1,sqr2,cosp));
                         _mm512_storeu_ps(&Tini[0] ,_mm512_fmadd_ps(sqi1,sqi2,cosp));
                 }


                    /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-74
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4174_zmm16r4(const __m512 mur,
                                          const __m512 mui,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 psi,
                                          __m512 * __restrict Toutr,
                                          __m512 * __restrict Touti) {

                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 divr,divi,sqr1,sqi1;
                         register __m512 mulr,muli,sqr2,sqi2;
                         register __m512 cosp,sinp,sin2p,t0r,t0i;
                         register __m512 _2sqr1,_2sqi1,t1r,t1i;
                         register __m512 numr,numi,denr,deni;
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&div);
                         cosp = xcosf(psi);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mulr,&muli);
                         sinp = xsinf(psi);
                         csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                         sin2p= _mm512_mul_ps(sinp,sinp);
                         //_2sqr1 = _mm512_mul_ps(_2,sqr1);
                         t0r    = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         //_2sqi1 = _mm512_mul_ps(_2,sqi1);
                         t0i    = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                         cmul_zmm16r4(sqr1,sqi1,t0r,t0i,&t1r,&t1i);
                         numr = _mm512_mul_ps(_2,t1r);
                         denr = _mm512_add_ps(cosp,t1r);
                         numi = _mm512_mul_ps(_2,t1i);
                         deni = _mm512_add_ps(cosp,t1i);
                         cdiv_zmm16r4(numr,numi,denr,deni,*Toutr,*Touti);
                 }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4174_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pmur,
                                            const float * __restrict __ATTR_ALIGN__(64) pmui,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                            float * __restrict __ATTR_ALIGN__(64) Toutr,
                                            float * __restrict __ATTR_ALIGN__(64) Touti) {

                         register __m512 mur  = _mm512_load_ps(&pmur[0]);
                         register __m512 mui  = _mm512_load_ps(&pmui[0]);
                         register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 divr,divi,sqr1,sqi1;
                         register __m512 mulr,muli,sqr2,sqi2;
                         register __m512 cosp,sinp,sin2p,t0r,t0i;
                         register __m512 _2sqr1,_2sqi1,t1r,t1i;
                         register __m512 numr,numi,denr,deni,resr,resi;
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&div);
                         cosp = xcosf(psi);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mulr,&muli);
                         sinp = xsinf(psi);
                         csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                         sin2p= _mm512_mul_ps(sinp,sinp);
                         //_2sqr1 = _mm512_mul_ps(_2,sqr1);
                         t0r    = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         //_2sqi1 = _mm512_mul_ps(_2,sqi1);
                         t0i    = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                         cmul_zmm16r4(sqr1,sqi1,t0r,t0i,&t1r,&t1i);
                         numr = _mm512_mul_ps(_2,t1r);
                         denr = _mm512_add_ps(cosp,t1r);
                         numi = _mm512_mul_ps(_2,t1i);
                         deni = _mm512_add_ps(cosp,t1i);
                         cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                         _mm512_store_ps(&Toutr[0], resr);
                         _mm512_store_ps(&Touti[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4174_zmm16r4_u(const float * __restrict pmur,
                                            const float * __restrict  pmui,
                                            const float * __restrict  pepsr,
                                            const float * __restrict  pepsi,
                                            const float * __restrict  ppsi,
                                            float * __restrict  Toutr,
                                            float * __restrict  Touti) {

                         register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                         register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 divr,divi,sqr1,sqi1;
                         register __m512 mulr,muli,sqr2,sqi2;
                         register __m512 cosp,sinp,sin2p,t0r,t0i;
                         register __m512 _2sqr1,_2sqi1,t1r,t1i;
                         register __m512 numr,numi,denr,deni,resr,resi;
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&div);
                         cosp = xcosf(psi);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mulr,&muli);
                         sinp = xsinf(psi);
                         csqrt_zmm16r4(divr,divi,&sqr1,&sqi1);
                         sin2p= _mm512_mul_ps(sinp,sinp);
                         //_2sqr1 = _mm512_mul_ps(_2,sqr1);
                         t0r    = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         //_2sqi1 = _mm512_mul_ps(_2,sqi1);
                         t0i    = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr2,&sqi2);
                         cmul_zmm16r4(sqr1,sqi1,t0r,t0i,&t1r,&t1i);
                         numr = _mm512_mul_ps(_2,t1r);
                         denr = _mm512_add_ps(cosp,t1r);
                         numi = _mm512_mul_ps(_2,t1i);
                         deni = _mm512_add_ps(cosp,t1i);
                         cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                         _mm512_storeu_ps(&Toutr[0], resr);
                         _mm512_storeu_ps(&Touti[0], resi);
                 }


                   /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-75
                       */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4175_zmm16r4(const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           const __m512 psi,
                                           __m512 * __restrict Toutr,
                                           __m512 * __restrict Touti) {

                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         register __m512 cosp,sinp,sin2p,sqr1,sqi1;
                         register __m512 sqr2,sqi2,_2cosp,divr,divi;
                         register __m512 mulr,muli,denr,deni,t0r,t0i;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_mul_ps(sinp,sinp);
                         _2cosp= _mm512_add_ps(cosp,cosp);
                         cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         t0i = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         cmul_zmm16r4(sqr1,sqi1,sqr2,sqi2,&denr,&deni);
                         denr = _mm512_add_ps(cosp,denr);
                         deni = _mm512_add_ps(cosp,deni);
                         *Toutr = _mm512_div_ps(_2cosp,denr);
                         *Touti = _mm512_div_ps(_2cosp,deni);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4175_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pmur,
                                            const float * __restrict __ATTR_ALIGN__(64) pmui,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                            float * __restrict __ATTR_ALIGN__(64) Toutr,
                                            float * __restrict __ATTR_ALIGN__(64) Touti) {

                         register __m512 mur  = _mm512_load_ps(&pmur[0]);
                         register __m512 mui  = _mm512_load_ps(&pmui[0]);
                         register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         register __m512 cosp,sinp,sin2p,sqr1,sqi1;
                         register __m512 sqr2,sqi2,_2cosp,divr,divi;
                         register __m512 mulr,muli,denr,deni,t0r,t0i;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_mul_ps(sinp,sinp);
                         _2cosp= _mm512_add_ps(cosp,cosp);
                         cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         t0i = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         cmul_zmm16r4(sqr1,sqi1,sqr2,sqi2,&denr,&deni);
                         denr = _mm512_add_ps(cosp,denr);
                         deni = _mm512_add_ps(cosp,deni);
                         _mm512_store_ps(&Toutr[0] ,_mm512_div_ps(_2cosp,denr));
                         _mm512_store_ps(&Touti[0] ,_mm512_div_ps(_2cosp,deni));
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4175_zmm16r4_u(const float * __restrict  pmur,
                                            const float * __restrict  pmui,
                                            const float * __restrict  pepsr,
                                            const float * __restrict  pepsi,
                                            const float * __restrict  ppsi,
                                            float * __restrict  Toutr,
                                            float * __restrict  Touti) {

                         register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                         register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         register __m512 cosp,sinp,sin2p,sqr1,sqi1;
                         register __m512 sqr2,sqi2,_2cosp,divr,divi;
                         register __m512 mulr,muli,denr,deni,t0r,t0i;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_mul_ps(sinp,sinp);
                         _2cosp= _mm512_add_ps(cosp,cosp);
                         cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         t0i = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         cmul_zmm16r4(sqr1,sqi1,sqr2,sqi2,&denr,&deni);
                         denr = _mm512_add_ps(cosp,denr);
                         deni = _mm512_add_ps(cosp,deni);
                         _mm512_storeu_ps(&Toutr[0] ,_mm512_div_ps(_2cosp,denr));
                         _mm512_storeu_ps(&Touti[0] ,_mm512_div_ps(_2cosp,deni));
                 }


                   /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-76
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4176_zmm16r4( const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           const __m512 psi,
                                           __m512 * __restrict Rinr,
                                           __m512 * __restrict Rini) {

                         register __m512 cosp,sinp,sin2p,divr,divi;
                         register __m512 mulr,muli,denr,deni,numr,numi;
                         register __m512 sqr1,sqi1,sqr2,sqi2,t0r,t0i;
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_ps(sinp,sinp);
                         cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_div_ps(sin2p,mulr));
                         t0i = _mm512_sub_ps(_1,_mm512_div_ps(sin2p,muli));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_ps(cosp,sqr2);
                         sqi2 = _mm512_mul_ps(cosp,sqi2);
                         numr = _mm512_sub_ps(sqr2,sqr1);
                         denr = _mm512_add_ps(sqr2,sqr1);
                         numi = _mm512_sub_ps(sqi2,sqi1);
                         deni = _mm512_add_ps(sqi2,sqi1);
                         cdiv_zmm16r4(numr,numi,denr,deni,*Rinr,*Rini);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4176_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) pmur,
                                             const float * __restrict __ATTR_ALIGN__(64) pmui,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                             const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                             float * __restrict __ATTR_ALIGN__(64) Rinr,
                                             float * __restrict __ATTR_ALIGN__(64) Rini) {

                         register __m512 mur  = _mm512_load_ps(&pmur[0]);
                         register __m512 mui  = _mm512_load_ps(&pmui[0]);
                         register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                         register __m512 cosp,sinp,sin2p,divr,divi;
                         register __m512 mulr,muli,denr,deni,numr,numi;
                         register __m512 sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_ps(sinp,sinp);
                         cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_div_ps(sin2p,mulr));
                         t0i = _mm512_sub_ps(_1,_mm512_div_ps(sin2p,muli));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_ps(cosp,sqr2);
                         sqi2 = _mm512_mul_ps(cosp,sqi2);
                         numr = _mm512_sub_ps(sqr2,sqr1);
                         denr = _mm512_add_ps(sqr2,sqr1);
                         numi = _mm512_sub_ps(sqi2,sqi1);
                         deni = _mm512_add_ps(sqi2,sqi1);
                         cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                         _mm512_store_ps(&Rinr[0], resr);
                         _mm512_store_ps(&Rini[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4176_zmm16r4_u( const float * __restrict  pmur,
                                             const float * __restrict  pmui,
                                             const float * __restrict  pepsr,
                                             const float * __restrict  pepsi,
                                             const float * __restrict  ppsi,
                                             float * __restrict  Rinr,
                                             float * __restrict  Rini) {

                         register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                         register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 cosp,sinp,sin2p,divr,divi;
                         register __m512 mulr,muli,denr,deni,numr,numi;
                         register __m512 sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_ps(sinp,sinp);
                         cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_div_ps(sin2p,mulr));
                         t0i = _mm512_sub_ps(_1,_mm512_div_ps(sin2p,muli));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_ps(cosp,sqr2);
                         sqi2 = _mm512_mul_ps(cosp,sqi2);
                         numr = _mm512_sub_ps(sqr2,sqr1);
                         denr = _mm512_add_ps(sqr2,sqr1);
                         numi = _mm512_sub_ps(sqi2,sqi1);
                         deni = _mm512_add_ps(sqi2,sqi1);
                         cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                         _mm512_storeu_ps(&Rinr[0], resr);
                         _mm512_storeu_ps(&Rini[0], resi);
                 }


                    /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-77
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4177_zmm16r4( const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           const __m512 psi,
                                           __m512 * __restrict Rinr,
                                           __m512 * __restrict Rini) {

                         register __m512 cosp,sinp,sin2p,divr,divi;
                         register __m512 mulr,muli,denr,deni,numr,numi;
                         register __m512 sqr1,sqi1,sqr2,sqi2,t0r,t0i;
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_ps(sinp,sinp);
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&divi);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         t0i = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_ps(sqr2,cosp);
                         sqi2 = _mm512_mul_ps(sqi2,cosp);
                         numr = _mm512_sub_ps(sqr2,sqr1);
                         denr = _mm512_add_ps(sqr2,sqr1);
                         numi = _mm512_sub_ps(sqi2,sqi1);
                         deni = _mm512_add_ps(sqi2,sqi1);
                         cdiv_zmm16r4(numr,numi,denr,deni,*Rinr,*Rini);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4177_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) pmur,
                                             const float * __restrict __ATTR_ALIGN__(64) pmui,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                             const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                             float * __restrict __ATTR_ALIGN__(64) Rinr,
                                             float * __restrict __ATTR_ALIGN__(64) Rini ) {

                         register __m512 mur  = _mm512_load_ps(&pmur[0]);
                         register __m512 mui  = _mm512_load_ps(&pmui[0]);
                         register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                         register __m512 cosp,sinp,sin2p,divr,divi;
                         register __m512 mulr,muli,denr,deni,numr,numi;
                         register __m512 sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_ps(sinp,sinp);
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&divi);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         t0i = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_ps(sqr2,cosp);
                         sqi2 = _mm512_mul_ps(sqi2,cosp);
                         numr = _mm512_sub_ps(sqr2,sqr1);
                         denr = _mm512_add_ps(sqr2,sqr1);
                         numi = _mm512_sub_ps(sqi2,sqi1);
                         deni = _mm512_add_ps(sqi2,sqi1);
                         cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                         _mm512_store_ps(&Rinr[0], resr);
                         _mm512_store_ps(&Rini[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4177_zmm16r4_u( const float * __restrict  pmur,
                                             const float * __restrict  pmui,
                                             const float * __restrict  pepsr,
                                             const float * __restrict  pepsi,
                                             const float * __restrict  ppsi,
                                             float * __restrict  Rinr,
                                             float * __restrict  Rini ) {

                         register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                         register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 cosp,sinp,sin2p,divr,divi;
                         register __m512 mulr,muli,denr,deni,numr,numi;
                         register __m512 sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512 _1 = _mm512_set1_ps(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_ps(sinp,sinp);
                         cdiv_zmm16r4(epsr,epsi,mur,mui,&divr,&divi);
                         cmul_zmm16r4(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_ps(_1,_mm512_mul_ps(mulr,sin2p));
                         t0i = _mm512_sub_ps(_1,_mm512_mul_ps(muli,sin2p));
                         csqrt_zmm16r4(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm16r4(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_ps(sqr2,cosp);
                         sqi2 = _mm512_mul_ps(sqi2,cosp);
                         numr = _mm512_sub_ps(sqr2,sqr1);
                         denr = _mm512_add_ps(sqr2,sqr1);
                         numi = _mm512_sub_ps(sqi2,sqi1);
                         deni = _mm512_add_ps(sqi2,sqi1);
                         cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                         _mm512_storeu_ps(&Rinr[0], resr);
                         _mm512_storeu_ps(&Rini[0], resi);
                }


                  /*
                          Specular rays reflection
                          Formula 4.1-64
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rext_f4164_zmm16r4(const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           __m512 * __restrict Rexr,
                                           __m512 * __restrict Rexi) {

                        register __m512 sqr1,sqi1,sqr2,sqi2;
                        register __m512 difr,difi,sumr,sumi;
                        csqrt_zmm16r4(mur,mui,&sqr1,sqi1);
                        csqrt_zmm16r4(epsr,epsi,&sqr2,&sqi2);
                        difr = _mm512_sub_ps(sqr1,sqr2);
                        sumr = _mm512_add_ps(sqr1,sqr2);
                        difi = _mm512_sub_ps(sqi1,sqi2);
                        sumi = _mm512_add_ps(sqi1,sqi2);
                        cdiv_zmm16r4(difr,difi,sumr,sumi,*Rexr,*Rexi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rext_f4164_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pmur,
                                             const float * __restrict __ATTR_ALIGN__(64) pmui,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                             float * __restrict __ATTR_ALIGN__(64) Rexr,
                                             float * __restrict __ATTR_ALIGN__(64) Rexi) {

                        register __m512 mur  = _mm512_load_ps(&pmur[0]);
                        register __m512 mui  = _mm512_load_ps(&pmui[0]);
                        register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512 sqr1,sqi1,sqr2,sqi2;
                        register __m512 difr,difi,sumr,sumi,resr,resi;
                        csqrt_zmm16r4(mur,mui,&sqr1,sqi1);
                        csqrt_zmm16r4(epsr,epsi,&sqr2,&sqi2);
                        difr = _mm512_sub_ps(sqr1,sqr2);
                        sumr = _mm512_add_ps(sqr1,sqr2);
                        difi = _mm512_sub_ps(sqi1,sqi2);
                        sumi = _mm512_add_ps(sqi1,sqi2);
                        cdiv_zmm16r4(difr,difi,sumr,sumi,&resr,&resi);
                        _mm512_store_ps(&Rexr[0], resr);
                        _mm512_store_ps(&Rexi[0], resi);
                }

                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rext_f4164_zmm16r4_u(const float * __restrict  pmur,
                                             const float * __restrict  pmui,
                                             const float * __restrict  pepsr,
                                             const float * __restrict  pepsi,
                                             float * __restrict  Rexr,
                                             float * __restrict  Rexi) {

                        register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                        register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 sqr1,sqi1,sqr2,sqi2;
                        register __m512 difr,difi,sumr,sumi,resr,resi;
                        csqrt_zmm16r4(mur,mui,&sqr1,sqi1);
                        csqrt_zmm16r4(epsr,epsi,&sqr2,&sqi2);
                        difr = _mm512_sub_ps(sqr1,sqr2);
                        sumr = _mm512_add_ps(sqr1,sqr2);
                        difi = _mm512_sub_ps(sqi1,sqi2);
                        sumi = _mm512_add_ps(sqi1,sqi2);
                        cdiv_zmm16r4(difr,difi,sumr,sumi,&resr,&resi);
                        _mm512_storeu_ps(&Rexr[0], resr);
                        _mm512_storeu_ps(&Rexi[0], resi);
                }


                  /*

                         Axial rays, when phi = 0
                         Formula 4.1-67
                    */

                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4167_zmm16r4( const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           __m512 * __restrict Tinr,
                                           __m512 * __restrict Tini) {

                          register __m512 sqr1,sqi1,sqr2,sqi2;
                          register __m512 sumr,sumi,mu2r,mu2i;
                          csqrt_zmm16r4(mur,mui,&sqr1,&sqi1);
                          mu2r = _mm512_add_ps(sqr1,sqr1);
                          mu2i = _mm512_add_ps(sqi1,sqi1);
                          csqrt_zmm16r4(epsr,epsi,&sqr2,&sqi2);
                          sumr = _mm512_add_ps(sqr1,sqr2);
                          sumi = _mm512_add_ps(sqi1,sqi2);
                          cdiv_zmm16r4(mu2r,mu2i,sumr,sumi,*Tinr,*Tini);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4167_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) pmur,
                                             const float * __restrict __ATTR_ALIGN__(64) pmui,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                             float * __restrict __ATTR_ALIGN__(64) Tinr,
                                             float * __restrict __ATTR_ALIGN__(64) Tini ) {

                          register __m512 mur  = _mm512_load_ps(&pmur[0]);
                          register __m512 mui  = _mm512_load_ps(&pmui[0]);
                          register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                          register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                          register __m512 sqr1,sqi1,sqr2,sqi2,resr,resi;
                          register __m512 sumr,sumi,mu2r,mu2i;
                          csqrt_zmm16r4(mur,mui,&sqr1,&sqi1);
                          mu2r = _mm512_add_ps(sqr1,sqr1);
                          mu2i = _mm512_add_ps(sqi1,sqi1);
                          csqrt_zmm16r4(epsr,epsi,&sqr2,&sqi2);
                          sumr = _mm512_add_ps(sqr1,sqr2);
                          sumi = _mm512_add_ps(sqi1,sqi2);
                          cdiv_zmm16r4(mu2r,mu2i,sumr,sumi,&resr,&resi);
                          _mm512_store_ps(&Tinr[0], resr);
                          _mm512_store_ps(&Tini[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4167_zmm16r4_u( const float * __restrict  pmur,
                                             const float * __restrict  pmui,
                                             const float * __restrict  pepsr,
                                             const float * __restrict  pepsi,
                                             float * __restrict  Tinr,
                                             float * __restrict  Tini ) {

                          register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                          register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                          register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                          register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                          register __m512 sqr1,sqi1,sqr2,sqi2,resr,resi;
                          register __m512 sumr,sumi,mu2r,mu2i;
                          csqrt_zmm16r4(mur,mui,&sqr1,&sqi1);
                          mu2r = _mm512_add_ps(sqr1,sqr1);
                          mu2i = _mm512_add_ps(sqi1,sqi1);
                          csqrt_zmm16r4(epsr,epsi,&sqr2,&sqi2);
                          sumr = _mm512_add_ps(sqr1,sqr2);
                          sumi = _mm512_add_ps(sqi1,sqi2);
                          cdiv_zmm16r4(mu2r,mu2i,sumr,sumi,&resr,&resi);
                          _mm512_storeu_ps(&Tinr[0], resr);
                          _mm512_storeu_ps(&Tini[0], resi);
                }


                  /*
                          Axial rays, when phi = 0
                          Formula 4.1-68
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4168_zmm16r4( const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           __m512 * __restrict Toutr,
                                           __m512 * __restrict Touti) {

                          register __m512 sqr1,sqi1,sqr2,sqi2;
                          register __m512 sumr,sumi,eps2r,eps2i;
                          csqrt_zmm16r4(epsr,epsi,&sqr1,&sqi1);
                          eps2r = _mm512_add_ps(sqr1,sqr1);
                          eps2i = _mm512_add_ps(sqi1,sqi1);
                          csqrt_zmm16r4(mur,mui,&sqr2,&sqi2);
                          sumr = _mm512_add_ps(sqr1,sqr2);
                          sumi = _mm512_add_ps(sqi1,sqi2);
                          cdiv_zmm16r4(eps2r,eps2i,sumr,sumi,*Toutr,*Touti);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4168_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) pmur,
                                             const float * __restrict __ATTR_ALIGN__(64) pmui,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                             float * __restrict __ATTR_ALIGN__(64) Toutr,
                                             float * __restrict __ATTR_ALIGN__(64) Touti ) {

                          register __m512 mur  = _mm512_load_ps(&pmur[0]);
                          register __m512 mui  = _mm512_load_ps(&pmui[0]);
                          register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                          register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                          register __m512 sqr1,sqi1,sqr2,sqi2;
                          register __m512 sumr,sumi,eps2r,eps2i,resr,resi;
                          csqrt_zmm16r4(epsr,epsi,&sqr1,&sqi1);
                          eps2r = _mm512_add_ps(sqr1,sqr1);
                          eps2i = _mm512_add_ps(sqi1,sqi1);
                          csqrt_zmm16r4(mur,mui,&sqr2,&sqi2);
                          sumr = _mm512_add_ps(sqr1,sqr2);
                          sumi = _mm512_add_ps(sqi1,sqi2);
                          cdiv_zmm16r4(eps2r,eps2i,sumr,sumi,&resr,&resi);
                          _mm512_store_ps(&Toutr[0], resr);
                          _mm512_store_ps(&Touti[0], resi);
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4168_zmm16r4_u(  const float * __restrict  pmur,
                                             const float * __restrict  pmui,
                                             const float * __restrict  pepsr,
                                             const float * __restrict  pepsi,
                                             float * __restrict  Toutr,
                                             float * __restrict Touti ) {

                          register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                          register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                          register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                          register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                          register __m512 sqr1,sqi1,sqr2,sqi2;
                          register __m512 sumr,sumi,eps2r,eps2i,resr,resi;
                          csqrt_zmm16r4(epsr,epsi,&sqr1,&sqi1);
                          eps2r = _mm512_add_ps(sqr1,sqr1);
                          eps2i = _mm512_add_ps(sqi1,sqi1);
                          csqrt_zmm16r4(mur,mui,&sqr2,&sqi2);
                          sumr = _mm512_add_ps(sqr1,sqr2);
                          sumi = _mm512_add_ps(sqi1,sqi2);
                          cdiv_zmm16r4(eps2r,eps2i,sumr,sumi,&resr,&resi);
                          _mm512_storeu_ps(&Toutr[0], resr);
                          _mm512_storeu_ps(&Touti[0], resi);
                }


                  /*
                          Axial rays, when phi = 0
                          Formula 4.1-69
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rint_f4169_zmm16r4(const __m512 mur,
                                           const __m512 mui,
                                           const __m512 epsr,
                                           const __m512 epsi,
                                           __m512 * __restrict Rintr,
                                           __m512 * __restrict Rinti) {
                        
                        register __m512 t0r,t0i;
                        const __m512 n1 = _mm512_mul_ps(-1.0f);
                        Rext_f4164_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                        *Rintr = _mm512_mul_ps(n1,t0r);
                        *Rinti = _mm512_mul_ps(n1,t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rint_f4169_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) pmur,
                                             const float * __restrict __ATTR_ALIGN__(64) pmui,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                             float * __restrict __ATTR_ALIGN__(64) Rintr,
                                             float * __restrict __ATTR_ALIGN__(64) Rinti) {
                        
                        register __m512 mur  = _mm512_load_ps(&pmur[0]);
                        register __m512 mui  = _mm512_load_ps(&pmui[0]);
                        register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512 t0r,t0i;
                        const __m512 n1 = _mm512_mul_ps(-1.0f);
                        Rext_f4164_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                        _mm512_store_ps(&Rintr[0] ,_mm512_mul_ps(n1,t0r));
                        _mm512_store_ps(&Rinti[0] ,_mm512_mul_ps(n1,t0i));
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rint_f4169_zmm16r4_u( const float * __restrict pmur,
                                              const float * __restrict  pmui,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              float * __restrict  Rintr,
                                              float * __restrict  Rinti) {
                        
                        register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                        register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 t0r,t0i;
                        const __m512 n1 = _mm512_mul_ps(-1.0f);
                        Rext_f4164_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                        _mm512_storeu_ps(&Rintr[0] ,_mm512_mul_ps(n1,t0r));
                        _mm512_storeu_ps(&Rinti[0] ,_mm512_mul_ps(n1,t0i));
                 }


                   /*
                       Backscatter widths in high-frequency limit.
                       Phi = 0, formula 4.1-91,for k1a>5.
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4191_zmm16r4(const __m512 a,
                                            const __m512 mur,
                                            const __m512 mui,
                                            const __m512 epsr,
                                            const __m512 epsi ) {

                          register __m512 t0r,t0i;
                          register __m512 cabs,rcs;
                          Rext_f4164_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                          cabs = cabs_zmm16r4(t0r,t0i);
                          rcs  = _mm512_mul_ps(cabs,_mm512_mul_ps(PI,a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4191_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi, ) {

                          register __m512 mur  = _mm512_load_ps(&pmur[0]);
                          register __m512 mui  = _mm512_load_ps(&pmui[0]);
                          register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                          register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                          register __m512 t0r,t0i;
                          register __m512 cabs,rcs;
                          Rext_f4164_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                          cabs = cabs_zmm16r4(t0r,t0i);
                          rcs  = _mm512_mul_ps(cabs,_mm512_mul_ps(PI,a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4191_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi, ) {

                          register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                          register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                          register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                          register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                          register __m512 t0r,t0i;
                          register __m512 cabs,rcs;
                          Rext_f4164_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                          cabs = cabs_zmm16r4(t0r,t0i);
                          rcs  = _mm512_mul_ps(cabs,_mm512_mul_ps(PI,a));
                          return (rcs);
                 }


                    /*
                         Bistatic scattering width (k0a0<<1, k1a0<<1), function of phi angle.
                         Formula 4.1-104
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41104_zmm16r4(const __m512 a0,
                                             const __m512 a1,
                                             const __m512 k0a0,
                                             const __m512 phi,
                                             const __m512 mu1r,
                                             const __m512 mu1i,
                                             const __m512 mu0r,
                                             const __m512 mu0i,
                                             const __m512 eps1r,
                                             const __m512 eps1i,
                                             const __m512 eps0r,
                                             const __m512 eps0i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 cosp,divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          cosp  = xcosf(phi);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                          div2i = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                          t1r   = _mm512_sub_ps(t0r,div2r);
                          t1i   = _mm512_sub_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41104_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa0,
                                               const float * __restrict __ATTR_ALIGN__(64) pa1,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                               const float * __restrict __ATTR_ALIGN__(64) pphi,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu1r,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu1i,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu0r,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu0i,
                                               const float * __restrict __ATTR_ALIGN__(64) peps1r,
                                               const float * __restrict __ATTR_ALIGN__(64) peps1i,
                                               const float * __restrict __ATTR_ALIGN__(64) peps0r,
                                               const float * __restrict __ATTR_ALIGN__(64) peps0i) {

                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 phi   = _mm512_load_ps(&pphi[0]);
                          register __m512 mu1r  = _mm512_load_ps(&pmu1r[0]);
                          register __m512 mu1i  = _mm512_load_ps(&pmu1i[0]);
                          register __m512 mu0r  = _mm512_load_ps(&pmu0r[0]);
                          register __m512 mu0i  = _mm512_load_ps(&pmu0i[0]);
                          register __m512 eps1r = _mm512_load_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_load_ps(&peps1i[0]);
                          register __m512 eps0r = _mm512_load_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_load_ps(&peps1r[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 cosp,divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          cosp  = xcosf(phi);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                          div2i = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                          t1r   = _mm512_sub_ps(t0r,div2r);
                          t1i   = _mm512_sub_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41104_zmm16r4_u(const float * __restrict  pa0,
                                               const float * __restrict  pa1,
                                               const float * __restrict  pk0a0,
                                               const float * __restrict  pphi,
                                               const float * __restrict  pmu1r,
                                               const float * __restrict  pmu1i,
                                               const float * __restrict  pmu0r,
                                               const float * __restrict  pmu0i,
                                               const float * __restrict  peps1r,
                                               const float * __restrict  peps1i,
                                               const float * __restrict  peps0r,
                                               const float * __restrict  peps0i) {

                          register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                          register __m512 a1    = _mm512_loadu_ps(&pa1[0]);
                          register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                          register __m512 phi   = _mm512_loadu_ps(&pphi[0]);
                          register __m512 mu1r  = _mm512_loadu_ps(&pmu1r[0]);
                          register __m512 mu1i  = _mm512_loadu_ps(&pmu1i[0]);
                          register __m512 mu0r  = _mm512_loadu_ps(&pmu0r[0]);
                          register __m512 mu0i  = _mm512_loadu_ps(&pmu0i[0]);
                          register __m512 eps1r = _mm512_loadu_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_loadu_ps(&peps1i[0]);
                          register __m512 eps0r = _mm512_loadu_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_loadu_ps(&peps1r[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 cosp,divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          cosp  = xcosf(phi);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                          div2i = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                          t1r   = _mm512_sub_ps(t0r,div2r);
                          t1i   = _mm512_sub_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                /*
                         Backscattering  width (k0a0<<1, k1a0<<1), phi = 0
                         Formula 4.1-105
                  */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41105_zmm16r4(const __m512 a0,
                                             const __m512 a1,
                                             const __m512 k0a0,
                                             const __m512 mu1r,
                                             const __m512 mu1i,
                                             const __m512 mu0r,
                                             const __m512 mu0i,
                                             const __m512 eps1r,
                                             const __m512 eps1i,
                                             const __m512 eps0r,
                                             const __m512 eps0i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,div2r);
                          div2i = _mm512_mul_ps(_2,div2r);
                          t1r   = _mm512_sub_ps(t0r,div2r);
                          t1i   = _mm512_sub_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41105_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa0,
                                               const float * __restrict __ATTR_ALIGN__(64) pa1,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu1r,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu1i,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu0r,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu0i,
                                               const float * __restrict __ATTR_ALIGN__(64) peps1r,
                                               const float * __restrict __ATTR_ALIGN__(64) peps1i,
                                               const float * __restrict __ATTR_ALIGN__(64) peps0r,
                                               const float * __restrict __ATTR_ALIGN__(64) peps0i) {

                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 mu1r  = _mm512_load_ps(&pmu1r[0]);
                          register __m512 mu1i  = _mm512_load_ps(&pmu1i[0]);
                          register __m512 mu0r  = _mm512_load_ps(&pmu0r[0]);
                          register __m512 mu0i  = _mm512_load_ps(&pmu0i[0]);
                          register __m512 eps1r = _mm512_load_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_load_ps(&peps1i[0]);
                          register __m512 eps0r = _mm512_load_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_load_ps(&peps1r[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,div2r);
                          div2i = _mm512_mul_ps(_2,div2r);
                          t1r   = _mm512_sub_ps(t0r,div2r);
                          t1i   = _mm512_sub_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41105_zmm16r4_u(const float * __restrict  pa0,
                                               const float * __restrict  pa1,
                                               const float * __restrict  pk0a0,
                                               const float * __restrict  pmu1r,
                                               const float * __restrict  pmu1i,
                                               const float * __restrict  pmu0r,
                                               const float * __restrict  pmu0i,
                                               const float * __restrict  peps1r,
                                               const float * __restrict  peps1i,
                                               const float * __restrict  peps0r,
                                               const float * __restrict  peps0i) {

                          register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                          register __m512 a1    = _mm512_loadu_ps(&pa1[0]);
                          register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                          register __m512 mu1r  = _mm512_loadu_ps(&pmu1r[0]);
                          register __m512 mu1i  = _mm512_loadu_ps(&pmu1i[0]);
                          register __m512 mu0r  = _mm512_loadu_ps(&pmu0r[0]);
                          register __m512 mu0i  = _mm512_loadu_ps(&pmu0i[0]);
                          register __m512 eps1r = _mm512_loadu_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_loadu_ps(&peps1i[0]);
                          register __m512 eps0r = _mm512_loadu_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_loadu_ps(&peps1r[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,div2r);
                          div2i = _mm512_mul_ps(_2,div2r);
                          t1r   = _mm512_sub_ps(t0r,div2r);
                          t1i   = _mm512_sub_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                /*
                      Forward scattering width (k0a0<<1, k1a0<<1), phi = pi.
                      Formula 4.1-106
                 */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41106_zmm16r4(const __m512 a0,
                                             const __m512 a1,
                                             const __m512 k0a0,
                                             const __m512 mu1r,
                                             const __m512 mu1i,
                                             const __m512 mu0r,
                                             const __m512 mu0i,
                                             const __m512 eps1r,
                                             const __m512 eps1i,
                                             const __m512 eps0r,
                                             const __m512 eps0i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,div2r);
                          div2i = _mm512_mul_ps(_2,div2r);
                          t1r   = _mm512_add_ps(t0r,div2r);
                          t1i   = _mm512_add_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41106_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa0,
                                               const float * __restrict __ATTR_ALIGN__(64) pa1,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu1r,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu1i,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu0r,
                                               const float * __restrict __ATTR_ALIGN__(64) pmu0i,
                                               const float * __restrict __ATTR_ALIGN__(64) peps1r,
                                               const float * __restrict __ATTR_ALIGN__(64) peps1i,
                                               const float * __restrict __ATTR_ALIGN__(64) peps0r,
                                               const float * __restrict __ATTR_ALIGN__(64) peps0i) {

                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 mu1r  = _mm512_load_ps(&pmu1r[0]);
                          register __m512 mu1i  = _mm512_load_ps(&pmu1i[0]);
                          register __m512 mu0r  = _mm512_load_ps(&pmu0r[0]);
                          register __m512 mu0i  = _mm512_load_ps(&pmu0i[0]);
                          register __m512 eps1r = _mm512_load_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_load_ps(&peps1i[0]);
                          register __m512 eps0r = _mm512_load_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_load_ps(&peps1r[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512 divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512 e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512 div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_ps(PI,a0);
                          k0a03 = _mm512_mul_ps(k0a0,
                                            _mm512_mul_ps(k0a0,k0a0));
                          frac  = _mm512_mul_ps(pia,_mm512_mul_ps(pi4,k0a03));
                          a1a0  = _mm512_div_ps(a1,a0);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_ps(_1,a1a0s);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          t0r   = _mm512_sub_ps(_mm512_mul_ps(divr,_1ma),_1);
                          t0i   = _mm512_sub_ps(_mm512_mul_ps(divi,_1ma),_1);
                          e1mr  = _mm512_mul_ps(eps1r,_1pa);
                          e0mr  = _mm512_mul_ps(eps0r,_1ma);
                          e1mi  = _mm512_mul_ps(eps1i,_1pa);
                          e0mi  = _mm512_mul_ps(eps0i,_1ma);
                          numr  = _mm512_sub_ps(e1mr,e0mr);
                          numi  = _mm512_sub_ps(e1mi,e0mi);
                          denr  = _mm512_add_ps(e1mr,e0mr);
                          deni  = _mm512_add_ps(e1mi,e0mi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_ps(_2,div2r);
                          div2i = _mm512_mul_ps(_2,div2r);
                          t1r   = _mm512_add_ps(t0r,div2r);
                          t1i   = _mm512_add_ps(t0i,div2i);
                          cabs  = cabs_zmm16r4(t1r,t1i);
                          rcs   = _mm512_mul_ps(frac,cabs);
                          return (rcs);
               }



                 /*
                       Hollow cylindrical shell.
                       Approximations for the low frequency region
                       (k0a0<<1, k1a0<<1).
                       Formula 4.1-124
                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41124_zmm16r4(const __m512 a1,
                                          const __m512 a0,
                                          const __m512 k0a0,
                                          const __m512 eps1r,
                                          const __m512 eps1i,
                                          const __m512 eps0r,
                                          const __m512 eps0i,
                                          __m512 * __restrict A0r,
                                          __m512 * __restrict A0i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512 t0r,t0i;
                          a1a0  = _mm512_div_ps(a1,a0);
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          fraci = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          cdiv_zmm16r4(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          divr = _mm512_sub_ps(divr,_1);
                          divi = _mm512_sub_ps(divi,_1);
                          t0r  = _mm512_mul_ps(divr,_1ma);
                          t0i  = _mm512_mul_ps(divi,_1ma);
                          cmul_zmm16r4(fracr,fraci,t0r,t0i,*A0r,*A0i);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41124_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  float * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps1r,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps1i,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps0r,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps0i,
                                            float * __restrict __ATTR_ALIGN__(64) A0r,
                                            float * __restrict __ATTR_ALIGN__(64) A0i) {

                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 eps1r = _mm512_load_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_load_ps(&peps1i[0]); 
                          register __m512 eps0r = _mm512_load_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_load_ps(&peps0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512 t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_ps(a1,a0);
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          fraci = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          cdiv_zmm16r4(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          divr = _mm512_sub_ps(divr,_1);
                          divi = _mm512_sub_ps(divi,_1);
                          t0r  = _mm512_mul_ps(divr,_1ma);
                          t0i  = _mm512_mul_ps(divi,_1ma);
                          cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                          _mm512_store_ps(&A0r[0], resr);
                          _mm512_store_ps(&A0i[0], resi);

               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41124_zmm16r4_u(const  float * __restrict  pa1,
                                          const  float * __restrict pa0,
                                          const  float * __restrict  pk0a0,
                                          const  float * __restrict  peps1r,
                                          const  float * __restrict  peps1i,
                                          const  float * __restrict  peps0r,
                                          const  float * __restrict  peps0i,
                                          float * __restrict  A0r,
                                          float * __restrict  A0i) {

                          register __m512 a1    = _mm512_loadu_ps(&pa1[0]);
                          register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                          register __m512 eps1r = _mm512_loadu_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_loadu_ps(&peps1i[0]); 
                          register __m512 eps0r = _mm512_loadu_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_loadu_ps(&peps0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512 t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_ps(a1,a0);
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          fraci = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          cdiv_zmm16r4(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          divr = _mm512_sub_ps(divr,_1);
                          divi = _mm512_sub_ps(divi,_1);
                          t0r  = _mm512_mul_ps(divr,_1ma);
                          t0i  = _mm512_mul_ps(divi,_1ma);
                          cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                          _mm512_storeu_ps(&A0r[0], resr);
                          _mm512_storeu_ps(&A0i[0], resi);

               }


                  /*

                       Hollow cylindrical shell.
                       Approximations for the low frequency region
                       (k0a0<<1, k1a0<<1).
                       Formula 4.1-126
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41126_zmm16r4(const __m512 a1,
                                          const __m512 a0,
                                          const __m512 k0a0,
                                          const __m512 mu1r,
                                          const __m512 mu1i,
                                          const __m512 mu0r,
                                          const __m512 mu0i,
                                          __m512 * __restrict B0r,
                                          __m512 * __restrict B0i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512 t0r,t0i;
                          a1a0  = _mm512_div_ps(a1,a0);
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          fraci = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          divr = _mm512_sub_ps(divr,_1);
                          divi = _mm512_sub_ps(divi,_1);
                          t0r  = _mm512_mul_ps(divr,_1ma);
                          t0i  = _mm512_mul_ps(divi,_1ma);
                          cmul_zmm16r4(fracr,fraci,t0r,t0i,*B0r,*B0i);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41126_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  float * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu1r,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu1i,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu0r,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu0i,
                                            float * __restrict __ATTR_ALIGN__(64) B0r,
                                            float * __restrict __ATTR_ALIGN__(64) B0i) {

                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 mu1r = _mm512_load_ps(&pmu1r[0]);
                          register __m512 mu1i = _mm512_load_ps(&pmu1i[0]); 
                          register __m512 mu0r = _mm512_load_ps(&pmu0r[0]);
                          register __m512 mu0i = _mm512_load_ps(&pmu0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512 t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_ps(a1,a0);
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          fraci = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          divr = _mm512_sub_ps(divr,_1);
                          divi = _mm512_sub_ps(divi,_1);
                          t0r  = _mm512_mul_ps(divr,_1ma);
                          t0i  = _mm512_mul_ps(divi,_1ma);
                          cmul_zmm16r4(fracr,fraci,t0r,t0i,*resr,*resi);
                          _mm512_store_ps(&B0r[0], resr);
                          _mm512_store_ps(&B0i[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41126_zmm16r4_u(const  float * __restrict  pa1,
                                            const  float * __restrict  pa0,
                                            const  float * __restrict  pk0a0,
                                            const  float * __restrict  mu1r,
                                            const  float * __restrict  mu1i,
                                            const  float * __restrict  mu0r,
                                            const  float * __restrict  mu0i,
                                            float * __restrict  B0r,
                                            float * __restrict  B0i) {

                          register __m512 a1    = _mm512_loadu_ps(&pa1[0]);
                          register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                          register __m512 mu1r = _mm512_loadu_ps(&pmu1r[0]);
                          register __m512 mu1i = _mm512_loadu_ps(&pmu1i[0]); 
                          register __m512 mu0r = _mm512_loadu_ps(&pmu0r[0]);
                          register __m512 mu0i = _mm512_loadu_ps(&pmu0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512 t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_ps(a1,a0);
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          fraci = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          divr = _mm512_sub_ps(divr,_1);
                          divi = _mm512_sub_ps(divi,_1);
                          t0r  = _mm512_mul_ps(divr,_1ma);
                          t0i  = _mm512_mul_ps(divi,_1ma);
                          cmul_zmm16r4(fracr,fraci,t0r,t0i,*resr,*resi);
                          _mm512_storeu_ps(&B0r[0], resr);
                          _mm512_storeu_ps(&B0i[0], resi);
               }


                  /*

                          Hollow cylindrical shell.
                          Approximations for the low frequency region
                          (k0a0<<1, k1a0<<1).
                           Formula 4.1-125
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41125_zmm16r4(const __m512 a1,
                                          const __m512 a0,
                                          const __m512 k0a0,
                                          const __m512 mu1r,
                                          const __m512 mu1i,
                                          const __m512 mu0r,
                                          const __m512 mu0i,
                                          __m512 * __restrict A1r,
                                          __m512 * __restrict A1i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512 divr,divi,divrs,divis,t1r,t1i;
                          register __m512 sqpr,sqpi,sqmr,sqmi,t0r,t0i;
                          register __m512 numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_ps(a1,a0);
                          fracr = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          cmul_zmm16r4(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_ps(divrs,_1);
                          divis = _mm512_sub_ps(divis,_1);
                          numr  = _mm512_mul_ps(divrs,_1ma);
                          numi  = _mm512_mul_ps(divis,_1ma);
                          t0r   = _mm512_add_ps(divr,_1);
                          t0i   = _mm512_add_ps(divi,_1);
                          cmul_zmm16r4(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_ps(divr,_1);
                          t1i   = _mm512_sub_ps(divi,_1);
                          cmul_zmm16r4(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_ps(sqmr,a1a02);
                          sqmi = _mm512_mul_ps(sqmi,a1a02);
                          denr = _mm512_sub_ps(sqpr,sqmr);
                          deni = _mm512_sub_ps(sqpi,sqmi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm16r4(fracr,fraci,ratr,rati,*A1r,*A1i);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41125_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  float * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu1r,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu1i,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu0r,
                                            const  float * __restrict __ATTR_ALIGN__(64) pmu0i,
                                            float * __restrict __ATTR_ALIGN__(64) A1r,
                                            float * __restrict __ATTR_ALIGN__(64) A1i) {

                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 mu1r = _mm512_load_ps(&pmu1r[0]);
                          register __m512 mu1i = _mm512_load_ps(&pmu1i[0]); 
                          register __m512 mu0r = _mm512_load_ps(&pmu0r[0]);
                          register __m512 mu0i = _mm512_load_ps(&pmu0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512 divr,divi,divrs,divis,t1r,t1i;
                          register __m512 sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512 numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_ps(a1,a0);
                          fracr = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          cmul_zmm16r4(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_ps(divrs,_1);
                          divis = _mm512_sub_ps(divis,_1);
                          numr  = _mm512_mul_ps(divrs,_1ma);
                          numi  = _mm512_mul_ps(divis,_1ma);
                          t0r   = _mm512_add_ps(divr,_1);
                          t0i   = _mm512_add_ps(divi,_1);
                          cmul_zmm16r4(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_ps(divr,_1);
                          t1i   = _mm512_sub_ps(divi,_1);
                          cmul_zmm16r4(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_ps(sqmr,a1a02);
                          sqmi = _mm512_mul_ps(sqmi,a1a02);
                          denr = _mm512_sub_ps(sqpr,sqmr);
                          deni = _mm512_sub_ps(sqpi,sqmi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm16r4(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_store_ps(&A1r[0], resr);
                          _mm512_store_ps(&A1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41125_zmm16r4_u(const  float * __restrict  pa1,
                                            const  float * __restrict  pa0,
                                            const  float * __restrict  pk0a0,
                                            const  float * __restrict  pmu1r,
                                            const  float * __restrict  pmu1i,
                                            const  float * __restrict  pmu0r,
                                            const  float * __restrict  pmu0i,
                                            float * __restrict  A1r,
                                            float * __restrict A1i) {

                          register __m512 a1    = _mm512_loadu_ps(&pa1[0]);
                          register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                          register __m512 mu1r = _mm512_loadu_ps(&pmu1r[0]);
                          register __m512 mu1i = _mm512_loadu_ps(&pmu1i[0]); 
                          register __m512 mu0r = _mm512_loadu_ps(&pmu0r[0]);
                          register __m512 mu0i = _mm512_loadu_ps(&pmu0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512 divr,divi,divrs,divis,t1r,t1i;
                          register __m512 sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512 numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_ps(a1,a0);
                          fracr = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          cdiv_zmm16r4(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          cmul_zmm16r4(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_ps(divrs,_1);
                          divis = _mm512_sub_ps(divis,_1);
                          numr  = _mm512_mul_ps(divrs,_1ma);
                          numi  = _mm512_mul_ps(divis,_1ma);
                          t0r   = _mm512_add_ps(divr,_1);
                          t0i   = _mm512_add_ps(divi,_1);
                          cmul_zmm16r4(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_ps(divr,_1);
                          t1i   = _mm512_sub_ps(divi,_1);
                          cmul_zmm16r4(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_ps(sqmr,a1a02);
                          sqmi = _mm512_mul_ps(sqmi,a1a02);
                          denr = _mm512_sub_ps(sqpr,sqmr);
                          deni = _mm512_sub_ps(sqpi,sqmi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm16r4(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_storeu_ps(&A1r[0], resr);
                          _mm512_storeu_ps(&A1i[0], resi);
                }


                 /*

                          Hollow cylindrical shell.
                          Approximations for the low frequency region
                          (k0a0<<1, k1a0<<1).
                           Formula 4.1-127
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41127_zmm16r4(const __m512 a1,
                                          const __m512 a0,
                                          const __m512 k0a0,
                                          const __m512 eps1r,
                                          const __m512 eps1i,
                                          const __m512 eps0r,
                                          const __m512 eps0i,
                                          __m512 * __restrict B1r,
                                          __m512 * __restrict B1i) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512 divr,divi,divrs,divis,t1r,t1i;
                          register __m512 sqpr,sqpi,sqmr,sqmi,t0r,t0i;
                          register __m512 numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_ps(a1,a0);
                          fracr = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          cdiv_zmm16r4(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          cmul_zmm16r4(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_ps(divrs,_1);
                          divis = _mm512_sub_ps(divis,_1);
                          numr  = _mm512_mul_ps(divrs,_1ma);
                          numi  = _mm512_mul_ps(divis,_1ma);
                          t0r   = _mm512_add_ps(divr,_1);
                          t0i   = _mm512_add_ps(divi,_1);
                          cmul_zmm16r4(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_ps(divr,_1);
                          t1i   = _mm512_sub_ps(divi,_1);
                          cmul_zmm16r4(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_ps(sqmr,a1a02);
                          sqmi = _mm512_mul_ps(sqmi,a1a02);
                          denr = _mm512_sub_ps(sqpr,sqmr);
                          deni = _mm512_sub_ps(sqpi,sqmi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm16r4(fracr,fraci,ratr,rati,*B1r,*B1i);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41127_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  float * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps1r,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps1i,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps0r,
                                            const  float * __restrict __ATTR_ALIGN__(64) peps0i,
                                            float * __restrict __ATTR_ALIGN__(64) B1r,
                                            float * __restrict __ATTR_ALIGN__(64) B1i) {

                          register __m512 a1    = _mm512_load_ps(&pa1[0]);
                          register __m512 a0    = _mm512_load_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                          register __m512 eps1r = _mm512_load_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_load_ps(&peps1i[0]); 
                          register __m512 eps0r = _mm512_load_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_load_ps(&peps0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512 divr,divi,divrs,divis,t1r,t1i;
                          register __m512 sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512 numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_ps(a1,a0);
                          fracr = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          cdiv_zmm16r4(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          cmul_zmm16r4(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_ps(divrs,_1);
                          divis = _mm512_sub_ps(divis,_1);
                          numr  = _mm512_mul_ps(divrs,_1ma);
                          numi  = _mm512_mul_ps(divis,_1ma);
                          t0r   = _mm512_add_ps(divr,_1);
                          t0i   = _mm512_add_ps(divi,_1);
                          cmul_zmm16r4(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_ps(divr,_1);
                          t1i   = _mm512_sub_ps(divi,_1);
                          cmul_zmm16r4(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_ps(sqmr,a1a02);
                          sqmi = _mm512_mul_ps(sqmi,a1a02);
                          denr = _mm512_sub_ps(sqpr,sqmr);
                          deni = _mm512_sub_ps(sqpi,sqmi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm16r4(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_store_ps(&B1r[0], resr);
                          _mm512_store_ps(&B1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41127_zmm16r4_u(const  float * __restrict  pa1,
                                            const  float * __restrict  pa0,
                                            const  float * __restrict  pk0a0,
                                            const  float * __restrict  peps1r,
                                            const  float * __restrict  peps1i,
                                            const  float * __restrict  peps0r,
                                            const  float * __restrict  peps0i,
                                            float * __restrict  B1r,
                                            float * __restrict  B1i) {

                          register __m512 a1    = _mm512_loadu_ps(&pa1[0]);
                          register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                          register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                          register __m512 eps1r = _mm512_loadu_ps(&peps1r[0]);
                          register __m512 eps1i = _mm512_loadu_ps(&peps1i[0]); 
                          register __m512 eps0r = _mm512_loadu_ps(&peps0r[0]);
                          register __m512 eps0i = _mm512_loadu_ps(&peps0i[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                          register __m512 k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512 divr,divi,divrs,divis,t1r,t1i;
                          register __m512 sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512 numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_ps(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_ps(a1,a0);
                          fracr = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a02));
                          a1a0s = _mm512_mul_ps(a1a0,a1a0);
                          _1ma  = _mm512_sub_ps(_1,a1a0s);
                          cdiv_zmm16r4(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          cmul_zmm16r4(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_ps(divrs,_1);
                          divis = _mm512_sub_ps(divis,_1);
                          numr  = _mm512_mul_ps(divrs,_1ma);
                          numi  = _mm512_mul_ps(divis,_1ma);
                          t0r   = _mm512_add_ps(divr,_1);
                          t0i   = _mm512_add_ps(divi,_1);
                          cmul_zmm16r4(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_ps(divr,_1);
                          t1i   = _mm512_sub_ps(divi,_1);
                          cmul_zmm16r4(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_ps(sqmr,a1a02);
                          sqmi = _mm512_mul_ps(sqmi,a1a02);
                          denr = _mm512_sub_ps(sqpr,sqmr);
                          deni = _mm512_sub_ps(sqpi,sqmi);
                          cdiv_zmm16r4(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm16r4(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_storeu_ps(&B1r[0], resr);
                          _mm512_storeu_ps(&B1i[0], resi);
                }


                    /*

                          Low-frequncy approximations (k0a<0.2)
                          Cylindrical Luneberg lens (k0a<0.2).
                          Formula 4.1-162
                     */

                    

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41162_zmm16r4(const __m512 k0a,
                                          __m512 * __restrict A0r,
                                          __m512 * __restrict A0i) {

                         const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         register __m512 k0a2,k0ah;
                         k0a2 = _mm512_mul_ps(k0a,k0a);
                         *A0r = Ir;
                         k0ah = _mm512_mul_ps(hlf,k0a2);
                         *A0i = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0ah));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41162_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            float * __restrict __ATTR_ALIGN__(64) A0r,
                                            float * __restrict __ATTR_ALIGN__(64) A0i) {

                   
                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         register __m512 k0a2,k0ah;
                         k0a2 = _mm512_mul_ps(k0a,k0a);
                         _mm512_store_ps(&A0r[0] ,Ir);
                         k0ah = _mm512_mul_ps(hlf,k0a2);
                         _mm512_store_ps(&A0i[0] ,_mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0ah)));
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41162_zmm16r4_u(const float * __restrict  pk0a,
                                            float * __restrict  A0r,
                                            float * __restrict  A0i) {

                   
                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         register __m512 k0a2,k0ah;
                         k0a2 = _mm512_mul_ps(k0a,k0a);
                         _mm512_storeu_ps(&A0r[0] ,Ir);
                         k0ah = _mm512_mul_ps(hlf,k0a2);
                         _mm512_storeu_ps(&A0i[0] ,_mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0ah)));
                }


                 

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41162_zmm16r4(__m512 * __restrict B0r,
                                          __m512 * __restrict B0i) {

                        *B0r = _mm512_setzero_ps();
                        *B0i = _mm512_setzero_ps();
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41162_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) B0r,
                                          float * __restrict __ATTR_ALIGN__(64) B0i) {

                        _mm512_store_ps(&B0r[0] ,_mm512_setzero_ps());
                        _mm512_store_ps(&B0i[0] ,_mm512_setzero_ps());
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41162_zmm16r4_u(float * __restrict  B0r,
                                          float * __restrict  B0i) {

                        _mm512_storeu_ps(&B0r[0] ,_mm512_setzero_ps());
                        _mm512_storeu_ps(&B0i[0] ,_mm512_setzero_ps());
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41162_zmm16r4(__m512 * __restrict A1r,
                                          __m512 * __restrict A1i) {

                        *A1r = _mm512_setzero_ps();
                        *A1i = _mm512_setzero_ps();
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41162_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) A1r,
                                          float * __restrict __ATTR_ALIGN__(64) A1i) {

                        _mm512_store_ps(&A1r[0] ,_mm512_setzero_ps());
                        _mm512_store_ps(&A1i[0] ,_mm512_setzero_ps());
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41162_zmm16r4_u(float * __restrict  A1r,
                                            float * __restrict  A1i) {

                        _mm512_storeu_ps(&A1r[0] ,_mm512_setzero_ps());
                        _mm512_storeu_ps(&A1i[0] ,_mm512_setzero_ps());
               }


                      __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41162_zmm16r4(const __m512 k0a,
                                          __m512 * __restrict B1r,
                                          __m512 * __restrict B1i) {

                         const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 c0 = _mm512_set1_ps(1.8992f);
                         register __m512 k0a2,k0ah;
                         k0a2 = _mm512_mul_ps(k0a,k0a);
                         *B1r = Ir;
                         k0ah = _mm512_mul_ps(c0,k0a2);
                         *B1i = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0ah));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41162_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            float * __restrict __ATTR_ALIGN__(64) B1r,
                                            float * __restrict __ATTR_ALIGN__(64) B1i) {

                   
                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 c0 = _mm512_set1_ps(1.8992f);
                         register __m512 k0a2,k0ah;
                         k0a2 = _mm512_mul_ps(k0a,k0a);
                         _mm512_store_ps(&B1r[0] ,Ir);
                         k0ah = _mm512_mul_ps(c0,k0a2);
                         _mm512_store_ps(&B1i[0] ,_mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0ah)));
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41162_zmm16r4_u(const float * __restrict  pk0a,
                                            float * __restrict  B1r,
                                            float * __restrict  B1i) {

                   
                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 c0 = _mm512_set1_ps(1.8992f);
                         register __m512 k0a2,k0ah;
                         k0a2 = _mm512_mul_ps(k0a,k0a);
                         _mm512_storeu_ps(&B1r[0] ,Ir);
                         k0ah = _mm512_mul_ps(c0,k0a2);
                         _mm512_storeu_ps(&B1i[0] ,_mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0ah)));
                }


                   /*
                          Low-frequncy approximations (k0a<0.2)
                          Cylindrical Luneberg lens (k0a<0.2).  
                          Scattering widths.
                          Formula 4.1-163
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41163_zmm16r4(const __m512 a,
                                             const __m512 k0a) {

                          const __m512 c0   = _mm512_set1_ps(0.0625f);
                          const __m512 pipi = _mm512_set1_ps( 9.869604401089358618834490999876f);
                          register __m512 rcs,t0,k0a3;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          t0   = _mm512_mul_ps(pipi,a);
                          rcs  = _mm512_mul_ps(k0a3,t0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41163_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]); 
                          const __m512 c0   = _mm512_set1_ps(0.0625f);
                          const __m512 pipi = _mm512_set1_ps( 9.869604401089358618834490999876f);
                          register __m512 rcs,t0,k0a3;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          t0   = _mm512_mul_ps(pipi,a);
                          rcs  = _mm512_mul_ps(k0a3,t0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41163_zmm16r4_u(const float * __restrict  pa,
                                             const float * __restrict  pk0a) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]); 
                          const __m512 c0   = _mm512_set1_ps(0.0625f);
                          const __m512 pipi = _mm512_set1_ps( 9.869604401089358618834490999876f);
                          register __m512 rcs,t0,k0a3;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          t0   = _mm512_mul_ps(pipi,a);
                          rcs  = _mm512_mul_ps(k0a3,t0);
                          return (rcs);
                 }


                    /*
                          Low-frequncy approximations (k0a<0.2)
                          Cylindrical Luneberg lens (k0a<0.2).  
                          Scattering widths.
                          Formula 4.1-164
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41164_zmm16r4(const __m512 a,
                                             const __m512 k0a,
                                             const __m512 phi) {

                          const __m512 c0   = _mm512_set1_ps(0.03607f);
                          const __m512 pipi = _mm512_set1_ps( 9.869604401089358618834490999876f);
                          register __m512 rcs,t0,cosp,k0a3;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(c0,_mm512_mul_ps(pipi,a));
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          rcs   = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,cos2p));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41164_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                               const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]); 
                          register __m512 phi = _mm512_load_ps(&pphi[0]);
                          const __m512 c0   = _mm512_set1_ps(0.03607f);
                          const __m512 pipi = _mm512_set1_ps( 9.869604401089358618834490999876f);
                          register __m512 rcs,t0,cosp,k0a3;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(c0,_mm512_mul_ps(pipi,a));
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          rcs   = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,cos2p));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f41164_zmm16r4_u(const float * __restrict  pa,
                                               const float * __restrict  pk0a,
                                               const float * __restrict  pphi) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]); 
                          register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                          const __m512 c0   = _mm512_set1_ps(0.03607f);
                          const __m512 pipi = _mm512_set1_ps( 9.869604401089358618834490999876f);
                          register __m512 rcs,t0,cosp,k0a3;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(c0,_mm512_mul_ps(pipi,a));
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          rcs   = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,cos2p));
                          return (rcs);
                 }


                  /*

                      Cylindrical Eaton-Lippman Lens, (k0a<0.2)
                      Formulae 4.1-165
                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41165_zmm16r4(const __m512 k0a,
                                          __m512 * __restrict A0r,
                                          __m512 * __restrict A0i) {

                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a2;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        *A0r = Ir;
                        *A0i = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a2));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41165_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            float * __restrict __ATTR_ALIGN__(64) A0r,
                                            float * __restrict __ATTR_ALIGN__(64) A0i) {

                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a2;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        _mm512_store_ps(&A0r[0], Ir);
                        _mm512_store_ps(&A0i[0], _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a2)));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41165_zmm16r4_u(const float * __restrict  pk0a,
                                            float * __restrict  A0r,
                                            float * __restrict  A0i) {

                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a2;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        _mm512_storeu_ps(&A0r[0], Ir);
                        _mm512_storeu_ps(&A0i[0], _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,k0a2)));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41165_zmm16r4(__m512 * __restrict A0r,
                                          __m512 * __restrict A0i) {

                        *A0r = Ir;
                        *A0i = Ir;
                }


                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41165_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) A0r,
                                            float * __restrict __ATTR_ALIGN__(64) A0i) {

                        _mm512_store_ps(&A0r[0], Ir);
                        _mm512_store_ps(&A0i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41165_zmm16r4_u(float * __restrict  A0r,
                                            float * __restrict  A0i) {

                        _mm512_storeu_ps(&A0r[0], Ir);
                        _mm512_storeu_ps(&A0i[0], Ir);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41165_zmm16r4(__m512 * __restrict B0r,
                                          __m512 * __restrict B0i) {

                        *B0r = Ir;
                        *B0i = Ir;
                }


                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41165_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) B0r,
                                            float * __restrict __ATTR_ALIGN__(64) B0i) {

                        _mm512_store_ps(&B0r[0], Ir);
                        _mm512_store_ps(&B0i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41165_zmm16r4_u(float * __restrict  B0r,
                                            float * __restrict  B0i) {

                        _mm512_storeu_ps(&B0r[0], Ir);
                        _mm512_storeu_ps(&B0i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41165_zmm16r4(const __m512 k0a,
                                          __m512 * __restrict B1r,
                                          __m512 * __restrict B1i) {

                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 c0  = _mm512_set1_ps(0.43616f);
                        register __m512 k0a2,t0;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        t0   = _mm512_mul_ps(c0,k0a2);
                        *B1r = Ir;
                        *B1i = _mm512_mul_ps(Ii,_mm512_mul_ps(pi4,t0));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41165_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                          float * __restrict __ATTR_ALIGN__(64) B1r,
                                          float * __restrict __ATTR_ALIGN__(64) B1i) {

                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 c0  = _mm512_set1_ps(0.43616f);
                        register __m512 k0a2,t0;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        t0   = _mm512_mul_ps(c0,k0a2);
                        _mm512_store_ps(&B1r[0] ,Ir);
                        _mm512_store_ps(&B1i[0] ,_mm512_mul_ps(Ii,_mm512_mul_ps(pi4,t0)));
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41165_zmm16r4_u(const float * __restrict  pk0a,
                                          float * __restrict  B1r,
                                          float * __restrict  B1i) {

                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 c0  = _mm512_set1_ps(0.43616f);
                        register __m512 k0a2,t0;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        t0   = _mm512_mul_ps(c0,k0a2);
                        _mm512_storeu_ps(&B1r[0] ,Ir);
                        _mm512_storeu_ps(&B1i[0] ,_mm512_mul_ps(Ii,_mm512_mul_ps(pi4,t0)));
                }


                  /*

                       Cylindrical Eaton-Lippman Lens, (k0a<0.2) 
                       Scattering widths.
                       Formula: 1.4-166
                   */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f14166_zmm16r4(const __m512 a,
                                             const __m512 k0a) {

                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 pip = _mm512_set1_ps(9.869604401089358618834490999876f);
                          register __m512 a4,k0a3,rcs;
                          a4   = _mm512_mul_ps(a,qtr);
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          rcs  = _mm512_mul_ps(k0a3,_mm512_mul_ps(pip,a4));
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f14166_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 pip = _mm512_set1_ps(9.869604401089358618834490999876f);
                          register __m512 a4,k0a3,rcs;
                          a4   = _mm512_mul_ps(a,qtr);
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          rcs  = _mm512_mul_ps(k0a3,_mm512_mul_ps(pip,a4));
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f14166_zmm16r4_u(const float * __restrict pa,
                                               const float * __restrict  pk0a) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 pip = _mm512_set1_ps(9.869604401089358618834490999876f);
                          register __m512 a4,k0a3,rcs;
                          a4   = _mm512_mul_ps(a,qtr);
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          rcs  = _mm512_mul_ps(k0a3,_mm512_mul_ps(pip,a4));
                          return (rcs);
               }


                  /*

                       Cylindrical Eaton-Lippman Lens, (k0a<0.2) 
                       Scattering widths.
                       Formula: 1.4-167
                   */
                 
 
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f14167_zmm16r4(const __m512 a,
                                             const __m512 k0a,
                                             const __m512 phi) {

                          const __m512 c0  = _mm512_set1_ps(0.19024f);
                          const __m512 pip = _mm512_set1_ps(9.869604401089358618834490999876f);
                          register __m512 cosp,cos2p,k0a3,rcs,t1;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          t1   = _mm512_mul_ps(_mm512_mul_ps(c0,pip),
                                               _mm512_mul_ps(a,k0a3));
                          cosp = xcosf(phi);
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          rcs   = _mm512_mul_ps(t1,cos2p);         
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f14167_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                               const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          register __m512 phi = _mm512_load_ps(&pphi[0]);
                          const __m512 c0  = _mm512_set1_ps(0.19024f);
                          const __m512 pip = _mm512_set1_ps(9.869604401089358618834490999876f);
                          register __m512 cosp,cos2p,k0a3,rcs,t1;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          t1   = _mm512_mul_ps(_mm512_mul_ps(c0,pip),
                                               _mm512_mul_ps(a,k0a3));
                          cosp = xcosf(phi);
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          rcs   = _mm512_mul_ps(t1,cos2p);         
                          return (rcs);
               }

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f14167_zmm16r4_u(const float * __restrict  pa,
                                               const float * __restrict  pk0a,
                                               const float * __restrict  pphi) {

                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                          const __m512 c0  = _mm512_set1_ps(0.19024f);
                          const __m512 pip = _mm512_set1_ps(9.869604401089358618834490999876f);
                          register __m512 cosp,cos2p,k0a3,rcs,t1;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          t1   = _mm512_mul_ps(_mm512_mul_ps(c0,pip),
                                               _mm512_mul_ps(a,k0a3));
                          cosp = xcosf(phi);
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          rcs   = _mm512_mul_ps(t1,cos2p);         
                          return (rcs);
               }


                  /*

                        Infinitely long cylinder.
                        Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                        TM-incident E-field.
                        Formula 4.2-48
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4248_zmm16r4(const __m512 E0r,
                                         const __m512 E0i,
                                         const __m512 psi,
                                         const __m512 phi,
                                         const __m512 k0,
                                         const __m512 z,
                                         const __m512 r,
                                         const __m512 a0,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui,
                                         __m512 * __restrict Ezr,
                                         __m512 * __restrict Ezi) {

                        const __m512 spi2 = _mm512_set1_ps(0.886226925452758013649083741671f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512 epsrp1,epsip1,epsrm1,epsim1;
                        register __m512 murp1,muip1,murm1,muim1,k0a02;
                        register __m512 mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512 ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512 frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_ps(k0,r);
                        k0z  = _mm512_mul_ps(k0,z);
                        k0a0 = _mm512_mul_ps(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_ps(hlf,_mm512_mul_ps(k0a0,k0a0));
                        epsrp1 = _mm512_add_ps(epsr,_1);
                        ear    = Ir;
                        epsip1 = _mm512_add_ps(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        scosps = _mm512_sqrt_ps(cosps);
                        epsim1 = _mm512_sub_ps(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        cos2ps= _mm512_mul_ps(cosps,cosps);
                        murm1  = _mm512_sub_ps(mur,_1);
                        cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_ps(mui,_1);
                        sin2ps= _mm512_mul_ps(sinps,sinps);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_ps(k0r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_ps(E0r,scosps);
                        fraci = _mm512_mul_ps(E0i,scosps);
                        cmul_zmm16r4(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_ps(spi2,_mm512_div_ps(frer,t1));
                        cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_ps(spi2,_mm512_div_ps(frei,t1));
                        cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_ps(epsrm1,cos2ps);
                        t0i = _mm512_mul_ps(epsim1,cos2ps);
                        numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                        cdiv_zmm16r4(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                        t1i = _mm512_mul_ps(_2,_mm512_mul_ps(div2i,cosp));
                        t2r = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0r,t1r));
                        t2i = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0i,t1i));
                        cmul_zmm16r4(div1r,div1i,t2r,t2i,*Ezr,*Ezi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4248_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pE0r,
                                           const float * __restrict __ATTR_ALIGN__(64) pE0i,
                                           const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pz,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pa0,
                                           const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                           const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                           const float * __restrict __ATTR_ALIGN__(64) pmur,
                                           const float * __restrict __ATTR_ALIGN__(64) pmui,
                                           float * __restrict __ATTR_ALIGN__(64) Ezr,
                                           float * __restrict __ATTR_ALIGN__(64) Ezi) {

                        register __m512 E0r = _mm512_load_ps(&pE0r[0]);
                        register __m512 E0i = _mm512_load_ps(&pE0i[0]);
                        register __m512 psi = _mm512_load_ps(&ppsi[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 z   = _mm512_load_ps(&pz[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 a0  = _mm512_load_ps(&pa0[0]);
                        register __m512 epsr= _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi= _mm512_load_ps(&pepsi[0]);
                        register __m512 pmur= _mm512_load_ps(&pmur[0]);
                        register __m512 pmui= _mm512_load_ps(&pmui[0]);
                        const __m512 spi2 = _mm512_set1_ps(0.886226925452758013649083741671f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512 epsrp1,epsip1,epsrm1,epsim1;
                        register __m512 murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512 mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512 ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512 frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_ps(k0,r);
                        k0z  = _mm512_mul_ps(k0,z);
                        k0a0 = _mm512_mul_ps(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_ps(hlf,_mm512_mul_ps(k0a0,k0a0));
                        epsrp1 = _mm512_add_ps(epsr,_1);
                        ear    = Ir;
                        epsip1 = _mm512_add_ps(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        scosps = _mm512_sqrt_ps(cosps);
                        epsim1 = _mm512_sub_ps(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        cos2ps= _mm512_mul_ps(cosps,cosps);
                        murm1  = _mm512_sub_ps(mur,_1);
                        cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_ps(mui,_1);
                        sin2ps= _mm512_mul_ps(sinps,sinps);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_ps(k0r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_ps(E0r,scosps);
                        fraci = _mm512_mul_ps(E0i,scosps);
                        cmul_zmm16r4(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_ps(spi2,_mm512_div_ps(frer,t1));
                        cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_ps(spi2,_mm512_div_ps(frei,t1));
                        cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_ps(epsrm1,cos2ps);
                        t0i = _mm512_mul_ps(epsim1,cos2ps);
                        numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                        cdiv_zmm16r4(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                        t1i = _mm512_mul_ps(_2,_mm512_mul_ps(div2i,cosp));
                        t2r = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0r,t1r));
                        t2i = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0i,t1i));
                        cmul_zmm16r4(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_store_ps(&Ezr[0], resr);
                        _mm512_store_ps(&Ezi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4248_zmm16r4_u(const float * __restrict  pE0r,
                                           const float * __restrict  pE0i,
                                           const float * __restrict  ppsi,
                                           const float * __restrict  pphi,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pz,
                                           const float * __restrict  pr,
                                           const float * __restrict  pa0,
                                           const float * __restrict  pepsr,
                                           const float * __restrict  pepsi,
                                           const float * __restrict  pmur,
                                           const float * __restrict  pmui,
                                           float * __restrict  Ezr,
                                           float * __restrict  Ezi) {

                        register __m512 E0r = _mm512_loadu_ps(&pE0r[0]);
                        register __m512 E0i = _mm512_loadu_ps(&pE0i[0]);
                        register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 z   = _mm512_loadu_ps(&pz[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 a0  = _mm512_loadu_ps(&pa0[0]);
                        register __m512 epsr= _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi= _mm512_loadu_ps(&pepsi[0]);
                        register __m512 pmur= _mm512_loadu_ps(&pmur[0]);
                        register __m512 pmui= _mm512_loadu_ps(&pmui[0]);
                        const __m512 spi2 = _mm512_set1_ps(0.886226925452758013649083741671f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512 epsrp1,epsip1,epsrm1,epsim1;
                        register __m512 murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512 mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512 ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512 frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_ps(k0,r);
                        k0z  = _mm512_mul_ps(k0,z);
                        k0a0 = _mm512_mul_ps(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_ps(hlf,_mm512_mul_ps(k0a0,k0a0));
                        epsrp1 = _mm512_add_ps(epsr,_1);
                        ear    = Ir;
                        epsip1 = _mm512_add_ps(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        scosps = _mm512_sqrt_ps(cosps);
                        epsim1 = _mm512_sub_ps(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        cos2ps= _mm512_mul_ps(cosps,cosps);
                        murm1  = _mm512_sub_ps(mur,_1);
                        cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_ps(mui,_1);
                        sin2ps= _mm512_mul_ps(sinps,sinps);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_ps(k0r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_ps(E0r,scosps);
                        fraci = _mm512_mul_ps(E0i,scosps);
                        cmul_zmm16r4(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_ps(spi2,_mm512_div_ps(frer,t1));
                        cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_ps(spi2,_mm512_div_ps(frei,t1));
                        cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_ps(epsrm1,cos2ps);
                        t0i = _mm512_mul_ps(epsim1,cos2ps);
                        numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                        cdiv_zmm16r4(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                        t1i = _mm512_mul_ps(_2,_mm512_mul_ps(div2i,cosp));
                        t2r = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0r,t1r));
                        t2i = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0i,t1i));
                        cmul_zmm16r4(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_storeu_ps(&Ezr[0], resr);
                        _mm512_storeu_ps(&Ezi[0], resi);
                }


               

                       /*

                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TM-incident H-field.
                         Formula 4.2-51
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hp_f4251_zmm16r4(const __m512 E0r,
                                         const __m512 E0i,
                                         const __m512 psi,
                                         const __m512 phi,
                                         const __m512 k0,
                                         const __m512 z,
                                         const __m512 r,
                                         const __m512 a0,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui,
                                         __m512 * __restrict Hpr,
                                         __m512 * __restrict Hpi) {

                        const __m512 e0u0 = _mm512_set1_ps(0.000007036193308495678572187302f);
                        const __m512 spi2 = _mm512_set1_ps(0.886226925452758013649083741671f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512 epsrp1,epsip1,epsrm1,epsim1,rat,spirat;
                        register __m512 murp1,muip1,murm1,muim1,k0a02;
                        register __m512 mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512 ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512 frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_ps(k0,r);
                        k0z  = _mm512_mul_ps(k0,z);
                        k0a0 = _mm512_mul_ps(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_ps(hlf,_mm512_mul_ps(k0a0,k0a0));
                        epsrp1 = _mm512_add_ps(epsr,_1);
                        ear    = Ir;
                        rat    = _mm512_sqrt_ps(e0u0);
                        epsip1 = _mm512_add_ps(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        scosps = _mm512_sqrt_ps(cosps);
                        epsim1 = _mm512_sub_ps(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        cos2ps= _mm512_mul_ps(cosps,cosps);
                        murm1  = _mm512_sub_ps(mur,_1);
                        cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_ps(mui,_1);
                        spirat = _mm512_mul_ps(spi2,rat);
                        sin2ps= _mm512_mul_ps(sinps,sinps);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_ps(E0r,scosps);
                        fraci = _mm512_mul_ps(E0i,scosps);
                        cmul_zmm16r4(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_ps(spirat,_mm512_div_ps(frer,t1));
                        cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_ps(spirat,_mm512_div_ps(frei,t1));
                        cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_ps(epsrm1,cos2ps);
                        t0i = _mm512_mul_ps(epsim1,cos2ps);
                        numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                        cdiv_zmm16r4(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                        t1i = _mm512_mul_ps(_2,_mm512_mul_ps(div2i,cosp));
                        t2r = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0r,t1r));
                        t2i = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0i,t1i));
                        div1r = _mm512_mul_ps(nIi,div1r);
                        div1i = _mm512_mul_ps(nIi,div1i);
                        cmul_zmm16r4(div1r,div1i,t2r,t2i,*Hpr,*Hpi);
                        
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hp_f4251_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) pE0r,
                                           const float * __restrict __ATTR_ALIGN__(64) pE0i,
                                           const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pz,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pa0,
                                           const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                           const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                           const float * __restrict __ATTR_ALIGN__(64) pmur,
                                           const float * __restrict __ATTR_ALIGN__(64) pmui,
                                           float * __restrict __ATTR_ALIGN__(64) Hpr,
                                           float * __restrict __ATTR_ALIGN__(64) Hpi) {


                        register __m512 E0r = _mm512_load_ps(&pE0r[0]);
                        register __m512 E0i = _mm512_load_ps(&pE0i[0]);
                        register __m512 psi = _mm512_load_ps(&ppsi[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 z   = _mm512_load_ps(&pz[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 a0  = _mm512_load_ps(&pa0[0]);
                        register __m512 epsr= _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi= _mm512_load_ps(&pepsi[0]);
                        register __m512 pmur= _mm512_load_ps(&pmur[0]);
                        register __m512 pmui= _mm512_load_ps(&pmui[0]);
                        const __m512 e0u0 = _mm512_set1_ps(0.000007036193308495678572187302f);
                        const __m512 spi2 = _mm512_set1_ps(0.886226925452758013649083741671f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512 epsrp1,epsip1,epsrm1,epsim1,rat,spirat;
                        register __m512 murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512 mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512 ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512 frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_ps(k0,r);
                        k0z  = _mm512_mul_ps(k0,z);
                        rat  = _mm512_div_ps(eps0,mu0);
                        k0a0 = _mm512_mul_ps(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_ps(hlf,_mm512_mul_ps(k0a0,k0a0));
                        epsrp1 = _mm512_add_ps(epsr,_1);
                        ear    = Ir;
                        rat    = _mm512_sqrt_ps(e0u0);
                        epsip1 = _mm512_add_ps(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        scosps = _mm512_sqrt_ps(cosps);
                        epsim1 = _mm512_sub_ps(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        cos2ps= _mm512_mul_ps(cosps,cosps);
                        murm1  = _mm512_sub_ps(mur,_1);
                        cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_ps(mui,_1);
                        spirat = _mm512_mul_ps(spi2,rat);
                        sin2ps= _mm512_mul_ps(sinps,sinps);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_ps(E0r,scosps);
                        fraci = _mm512_mul_ps(E0i,scosps);
                        cmul_zmm16r4(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_ps(spirat,_mm512_div_ps(frer,t1));
                        cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_ps(spirat,_mm512_div_ps(frei,t1));
                        cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_ps(epsrm1,cos2ps);
                        t0i = _mm512_mul_ps(epsim1,cos2ps);
                        numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                        cdiv_zmm16r4(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                        t1i = _mm512_mul_ps(_2,_mm512_mul_ps(div2i,cosp));
                        t2r = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0r,t1r));
                        t2i = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0i,t1i));
                        div1r = _mm512_mul_ps(nIi,div1r);
                        div1i = _mm512_mul_ps(nIi,div1i);
                        cmul_zmm16r4(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_store_ps(&Hpr[0], resr);
                        _mm512_store_ps(&Hpi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hp_f4251_zmm16r4_u(  const float * __restrict  pE0r,
                                           const float * __restrict  pE0i,
                                           const float * __restrict  ppsi,
                                           const float * __restrict pphi,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pz,
                                           const float * __restrict  pr,
                                           const float * __restrict  pa0,
                                           const float * __restrict  pepsr,
                                           const float * __restrict  pepsi,
                                           const float * __restrict  pmur,
                                           const float * __restrict  pmui,
                                           float * __restrict  Hpr,
                                           float * __restrict  Hpi) {


                        register __m512 E0r = _mm512_loadu_ps(&pE0r[0]);
                        register __m512 E0i = _mm512_loadu_ps(&pE0i[0]);
                        register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 z   = _mm512_loadu_ps(&pz[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 a0  = _mm512_loadu_ps(&pa0[0]);
                        register __m512 epsr= _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi= _mm512_loadu_ps(&pepsi[0]);
                        register __m512 pmur= _mm512_loadu_ps(&pmur[0]);
                        register __m512 pmui= _mm512_loadu_ps(&pmui[0]);
                        const __m512 e0u0 = _mm512_set1_ps(0.000007036193308495678572187302f);
                        const __m512 spi2 = _mm512_set1_ps(0.886226925452758013649083741671f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        register __m512 k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512 epsrp1,epsip1,epsrm1,epsim1,rat,spirat;
                        register __m512 murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512 mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512 ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512 frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_ps(k0,r);
                        k0z  = _mm512_mul_ps(k0,z);
                        rat  = _mm512_div_ps(eps0,mu0);
                        k0a0 = _mm512_mul_ps(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_ps(hlf,_mm512_mul_ps(k0a0,k0a0));
                        epsrp1 = _mm512_add_ps(epsr,_1);
                        ear    = Ir;
                        rat    = _mm512_sqrt_ps(e0u0);
                        epsip1 = _mm512_add_ps(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        scosps = _mm512_sqrt_ps(cosps);
                        epsim1 = _mm512_sub_ps(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        cos2ps= _mm512_mul_ps(cosps,cosps);
                        murm1  = _mm512_sub_ps(mur,_1);
                        cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_ps(mui,_1);
                        spirat = _mm512_mul_ps(spi2,rat);
                        sin2ps= _mm512_mul_ps(sinps,sinps);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_ps(E0r,scosps);
                        fraci = _mm512_mul_ps(E0i,scosps);
                        cmul_zmm16r4(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_ps(spirat,_mm512_div_ps(frer,t1));
                        cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_ps(spirat,_mm512_div_ps(frei,t1));
                        cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_ps(epsrm1,cos2ps);
                        t0i = _mm512_mul_ps(epsim1,cos2ps);
                        numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                        cdiv_zmm16r4(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_ps(_2,_mm512_mul_ps(div2r,cosp));
                        t1i = _mm512_mul_ps(_2,_mm512_mul_ps(div2i,cosp));
                        t2r = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0r,t1r));
                        t2i = _mm512_mul_ps(k0a02,_mm512_sub_ps(t0i,t1i));
                        div1r = _mm512_mul_ps(nIi,div1r);
                        div1i = _mm512_mul_ps(nIi,div1i);
                        cmul_zmm16r4(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_storeu_ps(&Hpr[0], resr);
                        _mm512_storeu_ps(&Hpi[0], resi);
                }


               
                     /*

                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TM-incident E-field.
                         Formula 4.2-49
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4249_zmm16r4(const __m512 E0r,
                                          const __m512 E0i,
                                          const __m512 k0z,
                                          const __m512 k0r,
                                          const __m512 k0a0,
                                          const __m512 psi,
                                          const __m512 phi,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 mur,
                                          const __m512 mui,
                                          __m512 * __restrict Ephr,
                                          __m512 * __restrict Ephi) {

                        const __m512 spi2 = _mm512_set1_ps(2.506628274631000502415765284811f);
                        const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i;
                        register __m512 cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512 emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_ps(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        t0    = _mm512_mul_ps(k0r,cosp);
                        ear   = Ir;
                        cmul_zmm16r4(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_ps(t0);
                        emum1r = _mm512_sub_ps(emum1r,_1);
                        emum1i = _mm512_sub_ps(emum1i,_1);
                        epsp1r = _mm512_add_ps(epsr,_1);
                        epsp1i = _mm512_add_ps(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosp,pi4));
                        sinpsp = _mm512_mul_ps(sinps,sinp);
                        eai    = t0;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_ps(fracr,den);
                        t0r   = _mm512_mul_ps(spi2,_mm512_mul_ps(fracr,k0a02));
                        fraci = _mm512_div_ps(fraci,den);
                        t0i   = _mm512_mul_ps(spi2,_mm512_mul_ps(fraci,k0a02));
                        cmul_zmm16r4(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm16r4(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_ps(sinpsp,t1r);
                        t1i = _mm512_mul_ps(sinpsp,t1i);
                        cmul_zmm16r4(fracr,fraci,t1r,t1i,*Ephr,*Ephi);
               }
                   


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4249_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pE0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pE0i,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui,
                                          float * __restrict __ATTR_ALIGN__(64) Ephr,
                                          float * __restrict __ATTR_ALIGN__(64) Ephi) {

                        register __m512 E0r   = _mm512_load_ps(&pE0r[0]);
                        register __m512 E0i   = _mm512_load_ps(&pE0i[0]);
                        register __m512 k0z   = _mm512_load_ps(&pk0z[0]);
                        register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                        register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                        register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                        register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_load_ps(&pmur[0]);
                        register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                        const __m512 spi2 = _mm512_set1_ps(2.506628274631000502415765284811f);
                        const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,resr,resi;
                        register __m512 cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512 emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_ps(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        t0    = _mm512_mul_ps(k0r,cosp);
                        ear   = Ir;
                        cmul_zmm16r4(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_ps(t0);
                        emum1r = _mm512_sub_ps(emum1r,_1);
                        emum1i = _mm512_sub_ps(emum1i,_1);
                        epsp1r = _mm512_add_ps(epsr,_1);
                        epsp1i = _mm512_add_ps(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosp,pi4));
                        sinpsp = _mm512_mul_ps(sinps,sinp);
                        eai    = t0;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_ps(fracr,den);
                        t0r   = _mm512_mul_ps(spi2,_mm512_mul_ps(fracr,k0a02));
                        fraci = _mm512_div_ps(fraci,den);
                        t0i   = _mm512_mul_ps(spi2,_mm512_mul_ps(fraci,k0a02));
                        cmul_zmm16r4(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm16r4(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_ps(sinpsp,t1r);
                        t1i = _mm512_mul_ps(sinpsp,t1i);
                        cmul_zmm16r4(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_store_ps(&Ephr[0],resr);
                        _mm512_store_ps(&Ephi[0],resi);
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4249_zmm16r4_u(const float * __restrict pE0r,
                                          const float * __restrict  pE0i,
                                          const float * __restrict pk0z,
                                          const float * __restrict  pk0r,
                                          const float * __restrict  pk0a0,
                                          const float * __restrict  ppsi,
                                          const float * __restrict  pphi,
                                          const float * __restrict  pepsr,
                                          const float * __restrict  pepsi,
                                          const float * __restrict  pmur,
                                          const float * __restrict  pmui,
                                          float * __restrict Ephr,
                                          float * __restrict  Ephi) {

                        register __m512 E0r   = _mm512_loadu_ps(&pE0r[0]);
                        register __m512 E0i   = _mm512_loadu_ps(&pE0i[0]);
                        register __m512 k0z   = _mm512_loadu_ps(&pk0z[0]);
                        register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                        register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                        register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                        register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                        const __m512 spi2 = _mm512_set1_ps(2.506628274631000502415765284811f);
                        const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,resr,resi;
                        register __m512 cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512 emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_ps(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        t0    = _mm512_mul_ps(k0r,cosp);
                        ear   = Ir;
                        cmul_zmm16r4(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_ps(t0);
                        emum1r = _mm512_sub_ps(emum1r,_1);
                        emum1i = _mm512_sub_ps(emum1i,_1);
                        epsp1r = _mm512_add_ps(epsr,_1);
                        epsp1i = _mm512_add_ps(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosp,pi4));
                        sinpsp = _mm512_mul_ps(sinps,sinp);
                        eai    = t0;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_ps(fracr,den);
                        t0r   = _mm512_mul_ps(spi2,_mm512_mul_ps(fracr,k0a02));
                        fraci = _mm512_div_ps(fraci,den);
                        t0i   = _mm512_mul_ps(spi2,_mm512_mul_ps(fraci,k0a02));
                        cmul_zmm16r4(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm16r4(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_ps(sinpsp,t1r);
                        t1i = _mm512_mul_ps(sinpsp,t1i);
                        cmul_zmm16r4(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_storeu_ps(&Ephr[0],resr);
                        _mm512_storeu_ps(&Ephi[0],resi);
               }



                 /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TM-incident H-field.
                         Formula 4.2-50

                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4250_zmm16r4( const __m512 E0r,
                                          const __m512 E0i,
                                          const __m512 k0z,
                                          const __m512 k0r,
                                          const __m512 k0a0,
                                          const __m512 psi,
                                          const __m512 phi,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 mur,
                                          const __m512 mui,
                                          __m512 * __restrict Hzr,
                                          __m512 * __restrict Hzi) {

                        const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                        const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,scosp;
                        register __m512 cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512 emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_ps(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        ear   = Ir;
                        cmul_zmm16r4(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_ps(k0r);
                        emum1r = _mm512_sub_ps(emum1r,_1);
                        emum1i = _mm512_sub_ps(emum1i,_1);
                        epsp1r = _mm512_add_ps(epsr,_1);
                        epsp1i = _mm512_add_ps(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosp,pi4));
                        scosp  = _mm512_sqrt_ps(cosp);
                        sinpsp = _mm512_mul_ps(sinps,sinp);
                        eai    = t0;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cer = _mm512_mul_ps(scosp,cer);
                        cei = _mm512_mul_ps(scosp,cei);
                        cmul_zmm16r4(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_ps(fracr,den);
                        t0r   = _mm512_mul_ps(e0u0,_mm512_mul_ps(fracr,k0a02));
                        fraci = _mm512_div_ps(fraci,den);
                        t0i   = _mm512_mul_ps(e0u0,_mm512_mul_ps(fraci,k0a02));
                        cmul_zmm16r4(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm16r4(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_ps(sinpsp,t1r);
                        t1i = _mm512_mul_ps(sinpsp,t1i);
                        cmul_zmm16r4(fracr,fraci,t1r,t1i,*Hzr,*Hzi);
               }


                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4250_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pE0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pE0i,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui,
                                          float * __restrict __ATTR_ALIGN__(64) Hzr,
                                          float * __restrict __ATTR_ALIGN__(64) Hzi ) {

                        register __m512 E0r   = _mm512_load_ps(&pE0r[0]);
                        register __m512 E0i   = _mm512_load_ps(&pE0i[0]);
                        register __m512 k0z   = _mm512_load_ps(&pk0z[0]);
                        register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                        register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                        register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                        register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_load_ps(&pmur[0]);
                        register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                        const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                        const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,scosp;
                        register __m512 cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli,resr,resi;
                        register __m512 emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_ps(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        ear   = Ir;
                        cmul_zmm16r4(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_ps(k0r);
                        emum1r = _mm512_sub_ps(emum1r,_1);
                        emum1i = _mm512_sub_ps(emum1i,_1);
                        epsp1r = _mm512_add_ps(epsr,_1);
                        epsp1i = _mm512_add_ps(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosp,pi4));
                        scosp  = _mm512_sqrt_ps(cosp);
                        sinpsp = _mm512_mul_ps(sinps,sinp);
                        eai    = t0;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cer = _mm512_mul_ps(scosp,cer);
                        cei = _mm512_mul_ps(scosp,cei);
                        cmul_zmm16r4(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_ps(fracr,den);
                        t0r   = _mm512_mul_ps(e0u0,_mm512_mul_ps(fracr,k0a02));
                        fraci = _mm512_div_ps(fraci,den);
                        t0i   = _mm512_mul_ps(e0u0,_mm512_mul_ps(fraci,k0a02));
                        cmul_zmm16r4(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm16r4(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_ps(sinpsp,t1r);
                        t1i = _mm512_mul_ps(sinpsp,t1i);
                        cmul_zmm16r4(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_store_ps(&Hzr[0], resr);
                        _mm512_store_ps(&Hzi[0], resi);
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4250_zmm16r4_u(const float * __restrict  pE0r,
                                          const float * __restrict pE0i,
                                          const float * __restrict  pk0z,
                                          const float * __restrict  pk0r,
                                          const float * __restrict  pk0a0,
                                          const float * __restrict  ppsi,
                                          const float * __restrict  pphi,
                                          const float * __restrict  pepsr,
                                          const float * __restrict  pepsi,
                                          const float * __restrict  pmur,
                                          const float * __restrict  pmui,
                                          float * __restrict  Hzr,
                                          float * __restrict  Hzi ) {

                        register __m512 E0r   = _mm512_loadu_ps(&pE0r[0]);
                        register __m512 E0i   = _mm512_loadu_ps(&pE0i[0]);
                        register __m512 k0z   = _mm512_loadu_ps(&pk0z[0]);
                        register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                        register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                        register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                        register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                        const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                        const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,scosp;
                        register __m512 cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli,resr,resi;
                        register __m512 emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_ps(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        ear   = Ir;
                        cmul_zmm16r4(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_ps(k0r);
                        emum1r = _mm512_sub_ps(emum1r,_1);
                        emum1i = _mm512_sub_ps(emum1i,_1);
                        epsp1r = _mm512_add_ps(epsr,_1);
                        epsp1i = _mm512_add_ps(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_ps(mur,_1);
                        muip1  = _mm512_add_ps(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosp,pi4));
                        scosp  = _mm512_sqrt_ps(cosp);
                        sinpsp = _mm512_mul_ps(sinps,sinp);
                        eai    = t0;
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cer = _mm512_mul_ps(scosp,cer);
                        cei = _mm512_mul_ps(scosp,cei);
                        cmul_zmm16r4(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_ps(fracr,den);
                        t0r   = _mm512_mul_ps(e0u0,_mm512_mul_ps(fracr,k0a02));
                        fraci = _mm512_div_ps(fraci,den);
                        t0i   = _mm512_mul_ps(e0u0,_mm512_mul_ps(fraci,k0a02));
                        cmul_zmm16r4(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm16r4(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_ps(sinpsp,t1r);
                        t1i = _mm512_mul_ps(sinpsp,t1i);
                        cmul_zmm16r4(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_storeu_ps(&Hzr[0], resr);
                        _mm512_storeu_ps(&Hzi[0], resi);
               }


                 /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident H-field.
                         Formula 4.2-52

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4252_zmm16r4(const __m512 H0r,
                                         const __m512 H0i,
                                         const __m512 psi,
                                         const __m512 phi,
                                         const __m512 k0r,
                                         const __m512 k0z,
                                         const __m512 k0a0,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui,
                                         __m512 * __restrict Hzr,
                                         __m512 * __restrict Hzi) {

                         const __m512 spi2 = _mm512_set1_ps(1.253314137315500251207882642406f);
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         register __m512 scosps,sinps,cosps,k0a02,sk0r;
                         register __m512 cos2ps,sin2ps,cosp,mul1r,mul1i;
                         register __m512 mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512 murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512 murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(hlf,_mm512_mul_ps(k0a,k0a));
                         scosps= _mm512_sqrt_ps(cosps);
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         t0    = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_ps(k0r);
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_ps(mur,_1);
                         fracr = _mm512_mul_ps(H0r,scosps);
                         sin2ps= _mm512_mul_ps(sinps,sinps);
                         muim1 = _mm512_sub_ps(mui,_1);
                         fraci = _mm512_mul_ps(H0i,scosps);
                         epsrm1= _mm512_sub_ps(epsr,_1);
                         epsim1= _mm512_sub_ps(epsi,_1);
                         murp1 = _mm512_add_ps(mur,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_ps(murm1,cos2ps);
                         t0r = _mm512_mul_ps(spi2,_mm512_div_ps(t0r,sk0r));
                         mucsi = _mm512_mul_ps(muim1,cos2ps);
                         t0i = _mm512_mul_ps(spi2,_mm512_div_ps(t0i,sk0r));
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_ps(t0r,k0a02);
                         t0i = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_ps(_2,_mm512_mul_ps(t1r,cosp));
                         t1i = _mm512_mul_ps(_2,_mm512_mul_ps(t1i,cosp));
                         mucsr = _mm512_sub_ps(mucsr,t1r);
                         mucsi = _mm512_sub_ps(mucsi,t1i);
                         cmul_zmm16r4(t0r,t0i,mucsr,mucsi,*Hzr,*Hzi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4252_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui,
                                          float * __restrict __ATTR_ALIGN__(64) Hzr,
                                          float * __restrict __ATTR_ALIGN__(64) Hzi) {

                         register __m512 H0r   = _mm512_load_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_load_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_load_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 spi2 = _mm512_set1_ps(1.253314137315500251207882642406f);
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         register __m512 scosps,sinps,cosps,k0a02,sk0r;
                         register __m512 cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512 mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512 murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512 murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(hlf,_mm512_mul_ps(k0a,k0a));
                         scosps= _mm512_sqrt_ps(cosps);
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         t0    = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_ps(k0r);
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_ps(mur,_1);
                         fracr = _mm512_mul_ps(H0r,scosps);
                         sin2ps= _mm512_mul_ps(sinps,sinps);
                         muim1 = _mm512_sub_ps(mui,_1);
                         fraci = _mm512_mul_ps(H0i,scosps);
                         epsrm1= _mm512_sub_ps(epsr,_1);
                         epsim1= _mm512_sub_ps(epsi,_1);
                         murp1 = _mm512_add_ps(mur,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_ps(murm1,cos2ps);
                         t0r = _mm512_mul_ps(spi2,_mm512_div_ps(t0r,sk0r));
                         mucsi = _mm512_mul_ps(muim1,cos2ps);
                         t0i = _mm512_mul_ps(spi2,_mm512_div_ps(t0i,sk0r));
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_ps(t0r,k0a02);
                         t0i = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_ps(_2,_mm512_mul_ps(t1r,cosp));
                         t1i = _mm512_mul_ps(_2,_mm512_mul_ps(t1i,cosp));
                         mucsr = _mm512_sub_ps(mucsr,t1r);
                         mucsi = _mm512_sub_ps(mucsi,t1i);
                         cmul_zmm16r4(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_store_ps(&Hzr[0], resr);
                         _mm512_store_ps(&Hzi[0], resi);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4252_zmm16r4_u(const float * __restrict  pH0r,
                                          const float * __restrict  pH0i,
                                          const float * __restrict  pk0z,
                                          const float * __restrict  pk0r,
                                          const float * __restrict  pk0a0,
                                          const float * __restrict  ppsi,
                                          const float * __restrict  pphi,
                                          const float * __restrict  pepsr,
                                          const float * __restrict  pepsi,
                                          const float * __restrict pmur,
                                          const float * __restrict  pmui,
                                          float * __restrict  Hzr,
                                          float * __restrict Hzi) {

                         register __m512 H0r   = _mm512_loadu_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_loadu_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_loadu_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 spi2 = _mm512_set1_ps(1.253314137315500251207882642406f);
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         register __m512 scosps,sinps,cosps,k0a02,sk0r;
                         register __m512 cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512 mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512 murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512 murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(hlf,_mm512_mul_ps(k0a,k0a));
                         scosps= _mm512_sqrt_ps(cosps);
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         t0    = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_ps(k0r);
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_ps(mur,_1);
                         fracr = _mm512_mul_ps(H0r,scosps);
                         sin2ps= _mm512_mul_ps(sinps,sinps);
                         muim1 = _mm512_sub_ps(mui,_1);
                         fraci = _mm512_mul_ps(H0i,scosps);
                         epsrm1= _mm512_sub_ps(epsr,_1);
                         epsim1= _mm512_sub_ps(epsi,_1);
                         murp1 = _mm512_add_ps(mur,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_ps(murm1,cos2ps);
                         t0r = _mm512_mul_ps(spi2,_mm512_div_ps(t0r,sk0r));
                         mucsi = _mm512_mul_ps(muim1,cos2ps);
                         t0i = _mm512_mul_ps(spi2,_mm512_div_ps(t0i,sk0r));
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_ps(t0r,k0a02);
                         t0i = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_ps(_2,_mm512_mul_ps(t1r,cosp));
                         t1i = _mm512_mul_ps(_2,_mm512_mul_ps(t1i,cosp));
                         mucsr = _mm512_sub_ps(mucsr,t1r);
                         mucsi = _mm512_sub_ps(mucsi,t1i);
                         cmul_zmm16r4(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_storeu_ps(&Hzr[0], resr);
                         _mm512_storeu_ps(&Hzi[0], resi);
                }


                    /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident E-field.
                         Formula 4.2-55

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4255_zmm16r4(const __m512 H0r,
                                         const __m512 H0i,
                                         const __m512 psi,
                                         const __m512 phi,
                                         const __m512 k0r,
                                         const __m512 k0z,
                                         const __m512 k0a0,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui,
                                         __m512 * __restrict Epr,
                                         __m512 * __restrict Epi) {

                         const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         register __m512 sinps,cosps,k0a02,sk0r;
                         register __m512 cos2ps,sin2ps,cosp,mul1r,mul1i;
                         register __m512 mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512 murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512 murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(hlf,_mm512_mul_ps(k0a,k0a));
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         t0    = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_ps(mur,_1);
                         fracr = H0r;
                         sin2ps= _mm512_mul_ps(sinps,sinps);
                         muim1 = _mm512_sub_ps(mui,_1);
                         fraci = H0i;
                         epsrm1= _mm512_sub_ps(epsr,_1);
                         epsim1= _mm512_sub_ps(epsi,_1);
                         murp1 = _mm512_add_ps(mur,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_ps(murm1,cos2ps);
                         t0r = _mm512_mul_ps(e0u0,_mm512_div_ps(t0r,sk0r));
                         mucsi = _mm512_mul_ps(muim1,cos2ps);
                         t0i = _mm512_mul_ps(e0u0,_mm512_div_ps(t0i,sk0r));
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_ps(t0r,k0a02);
                         t0i = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_ps(_2,_mm512_mul_ps(t1r,cosp));
                         t1i = _mm512_mul_ps(_2,_mm512_mul_ps(t1i,cosp));
                         mucsr = _mm512_sub_ps(mucsr,t1r);
                         mucsi = _mm512_sub_ps(mucsi,t1i);
                         cmul_zmm16r4(t0r,t0i,mucsr,mucsi,*Epr,*Epi);
                }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4255_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui,
                                          float * __restrict __ATTR_ALIGN__(64) Epr,
                                          float * __restrict __ATTR_ALIGN__(64) Epi) {

                         register __m512 H0r   = _mm512_load_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_load_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_load_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         register __m512 sinps,cosps,k0a02,sk0r;
                         register __m512 cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512 mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512 murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512 murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(hlf,_mm512_mul_ps(k0a,k0a));
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         t0    = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_ps(mur,_1);
                         fracr = H0r;
                         sin2ps= _mm512_mul_ps(sinps,sinps);
                         muim1 = _mm512_sub_ps(mui,_1);
                         fraci = H0i;
                         epsrm1= _mm512_sub_ps(epsr,_1);
                         epsim1= _mm512_sub_ps(epsi,_1);
                         murp1 = _mm512_add_ps(mur,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_ps(murm1,cos2ps);
                         t0r = _mm512_mul_ps(e0u0,_mm512_div_ps(t0r,sk0r));
                         mucsi = _mm512_mul_ps(muim1,cos2ps);
                         t0i = _mm512_mul_ps(e0u0,_mm512_div_ps(t0i,sk0r));
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_ps(t0r,k0a02);
                         t0i = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_ps(_2,_mm512_mul_ps(t1r,cosp));
                         t1i = _mm512_mul_ps(_2,_mm512_mul_ps(t1i,cosp));
                         mucsr = _mm512_sub_ps(mucsr,t1r);
                         mucsi = _mm512_sub_ps(mucsi,t1i);
                         cmul_zmm16r4(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_store_ps(&Epr[0], resr);
                         _mm512_store_ps(&Epi[0], resi);
                }


                      __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4255_zmm16r4_u(const float * __restrict pH0r,
                                          const float * __restrict  pH0i,
                                          const float * __restrict  pk0z,
                                          const float * __restrict  pk0r,
                                          const float * __restrict  pk0a0,
                                          const float * __restrict  ppsi,
                                          const float * __restrict  pphi,
                                          const float * __restrict  pepsr,
                                          const float * __restrict  pepsi,
                                          const float * __restrict  pmur,
                                          const float * __restrict  pmui,
                                          float * __restrict  Epr,
                                          float * __restrict  Epi) {

                         register __m512 H0r   = _mm512_loadu_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_loadu_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_loadu_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         register __m512 sinps,cosps,k0a02,sk0r;
                         register __m512 cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512 mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512 murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512 murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(hlf,_mm512_mul_ps(k0a,k0a));
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         t0    = _mm512_fmadd_ps(k0z,sinps,_mm512_fmadd_ps(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_ps(mur,_1);
                         fracr = H0r;
                         sin2ps= _mm512_mul_ps(sinps,sinps);
                         muim1 = _mm512_sub_ps(mui,_1);
                         fraci = H0i;
                         epsrm1= _mm512_sub_ps(epsr,_1);
                         epsim1= _mm512_sub_ps(epsi,_1);
                         murp1 = _mm512_add_ps(mur,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_ps(murm1,cos2ps);
                         t0r = _mm512_mul_ps(e0u0,_mm512_div_ps(t0r,sk0r));
                         mucsi = _mm512_mul_ps(muim1,cos2ps);
                         t0i = _mm512_mul_ps(e0u0,_mm512_div_ps(t0i,sk0r));
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_ps(t0r,k0a02);
                         t0i = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_ps(_2,_mm512_mul_ps(t1r,cosp));
                         t1i = _mm512_mul_ps(_2,_mm512_mul_ps(t1i,cosp));
                         mucsr = _mm512_sub_ps(mucsr,t1r);
                         mucsi = _mm512_sub_ps(mucsi,t1i);
                         cmul_zmm16r4(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_storeu_ps(&Epr[0], resr);
                         _mm512_storeu_ps(&Epi[0], resi);
                }


                    /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident H-field.
                         Formula 4.2-53

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hph_f4253_zmm16r4(const __m512 H0r,
                                          const __m512 H0i,
                                          const __m512 k0z,
                                          const __m512 k0r,
                                          const __m512 psi,
                                          const __m512 phi,
                                          const __m512 k0a0,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 mur,
                                          const __m512 mui,
                                          __m512 * __restrict Hpr,
                                          __m512 * __restrict Hpi) {

                         const __m512 s2pi = _mm512_set1_ps(2.506628274631000502415765284811f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 sinps,cosps,scpk0r,sinph;
                         register __m512 k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512 spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosps,pi4));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                         cmul_zmm16r4(H0r,H0i,cer,cei,&fracr,&fraci);
                         spsph = _mm512_mul_ps(sinps,sinph);
                         murp1 = _mm512_add_ps(mur,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         t0r   = _mm512_mul_ps(s2pi,_mm512_div_ps(fracr,scpk0r));
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         t0i   = _mm512_mul_ps(s2pi,_mm512_div_ps(fraci,scpk0r));
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,spsph);
                         t0r  = _mm512_mul_ps(t0r,k0a02);
                         divi = _mm512_mul_ps(divi,spsph);
                         t0i  = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(t0r,t0i,divr,divi,*Hpr,*Hpi);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hph_f4253_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui,
                                          float * __restrict __ATTR_ALIGN__(64) Hpr,
                                          float * __restrict __ATTR_ALIGN__(64) Hpi) {

                         register __m512 H0r   = _mm512_load_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_load_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_load_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 s2pi = _mm512_set1_ps(2.506628274631000502415765284811f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512 k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512 spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosps,pi4));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                         cmul_zmm16r4(H0r,H0i,cer,cei,&fracr,&fraci);
                         spsph = _mm512_mul_ps(sinps,sinph);
                         murp1 = _mm512_add_ps(mur,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         t0r   = _mm512_mul_ps(s2pi,_mm512_div_ps(fracr,scpk0r));
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         t0i   = _mm512_mul_ps(s2pi,_mm512_div_ps(fraci,scpk0r));
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,spsph);
                         t0r  = _mm512_mul_ps(t0r,k0a02);
                         divi = _mm512_mul_ps(divi,spsph);
                         t0i  = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_store_ps(&Hpr[0], resr);
                         _mm512_store_ps(&Hpi[0], resi);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hph_f4253_zmm16r4_u(const float * __restrict  pH0r,
                                          const float * __restrict pH0i,
                                          const float * __restrict  pk0z,
                                          const float * __restrict  pk0r,
                                          const float * __restrict  pk0a0,
                                          const float * __restrict  ppsi,
                                          const float * __restrict  pphi,
                                          const float * __restrict  pepsr,
                                          const float * __restrict  pepsi,
                                          const float * __restrict  pmur,
                                          const float * __restrict  pmui,
                                          float * __restrict  Hpr,
                                          float * __restrict  Hpi) {

                         register __m512 H0r   = _mm512_loadu_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_loadu_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_loadu_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 s2pi = _mm512_set1_ps(2.506628274631000502415765284811f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512 k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512 spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosps,pi4));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_ps(_mm512_mul_ps(k0r,cosps));
                         cmul_zmm16r4(H0r,H0i,cer,cei,&fracr,&fraci);
                         spsph = _mm512_mul_ps(sinps,sinph);
                         murp1 = _mm512_add_ps(mur,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         t0r   = _mm512_mul_ps(s2pi,_mm512_div_ps(fracr,scpk0r));
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         t0i   = _mm512_mul_ps(s2pi,_mm512_div_ps(fraci,scpk0r));
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,spsph);
                         t0r  = _mm512_mul_ps(t0r,k0a02);
                         divi = _mm512_mul_ps(divi,spsph);
                         t0i  = _mm512_mul_ps(t0i,k0a02);
                         cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_storeu_ps(&Hpr[0], resr);
                         _mm512_storeu_ps(&Hpi[0], resi);
                 }


                   /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident E-field.
                         Formula 4.2-54

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4254_zmm16r4( const __m512 H0r,
                                          const __m512 H0i,
                                          const __m512 k0z,
                                          const __m512 k0r,
                                          const __m512 psi,
                                          const __m512 phi,
                                          const __m512 k0a0,
                                          const __m512 epsr,
                                          const __m512 epsi,
                                          const __m512 mur,
                                          const __m512 mui,
                                          __m512 * __restrict Ezr,
                                          __m512 * __restrict Ezi) {
                         
                         const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 sinps,cosps,scpk0r,sinph;
                         register __m512 k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512 spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosps,pi4));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_ps(k0r);
                         cmul_zmm16r4(H0r,H0i,cer,cei,&fracr,&fraci);
                         fracr = _mm512_mul_ps(fracr,cosps);
                         spsph = _mm512_mul_ps(sinps,sinph);
                         fraci = _mm512_mul_ps(fraci,cosps);
                         murp1 = _mm512_add_ps(mur,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         t0r   = _mm512_mul_ps(e0u0,_mm512_div_ps(fracr,scpk0r));
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         t0i   = _mm512_mul_ps(e0u0,_mm512_div_ps(fraci,scpk0r));
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,spsph);
                         t0r  = _mm512_mul_ps(Ii,_mm512_mul_ps(t0r,k0a02));
                         divi = _mm512_mul_ps(divi,spsph);
                         t0i  = _mm512_mul_ps(Ii,_mm512_mul_ps(t0i,k0a02));
                         cmul_zmm16r4(t0r,t0i,divr,divi,*Ezr,*Ezi);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4254_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui,
                                          float * __restrict __ATTR_ALIGN__(64) Hzr,
                                          float * __restrict __ATTR_ALIGN__(64) Hzi ) {
                         
                         register __m512 H0r   = _mm512_load_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_load_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_load_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512 k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512 spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosps,pi4));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_ps(k0r);
                         cmul_zmm16r4(H0r,H0i,cer,cei,&fracr,&fraci);
                         fracr = _mm512_mul_ps(fracr,cosps);
                         spsph = _mm512_mul_ps(sinps,sinph);
                         fraci = _mm512_mul_ps(fraci,cosps);
                         murp1 = _mm512_add_ps(mur,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         t0r   = _mm512_mul_ps(e0u0,_mm512_div_ps(fracr,scpk0r));
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         t0i   = _mm512_mul_ps(e0u0,_mm512_div_ps(fraci,scpk0r));
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,spsph);
                         t0r  = _mm512_mul_ps(Ii,_mm512_mul_ps(t0r,k0a02));
                         divi = _mm512_mul_ps(divi,spsph);
                         t0i  = _mm512_mul_ps(Ii,_mm512_mul_ps(t0i,k0a02));
                         cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_store_ps(&Hzr[0], resr);
                         _mm512_store_ps(&Hzi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4254_zmm16r4_u(const float * __restrict  pH0r,
                                          const float * __restrict pH0i,
                                          const float * __restrict  pk0z,
                                          const float * __restrict  pk0r,
                                          const float * __restrict  pk0a0,
                                          const float * __restrict  ppsi,
                                          const float * __restrict  pphi,
                                          const float * __restrict  pepsr,
                                          const float * __restrict  pepsi,
                                          const float * __restrict  pmur,
                                          const float * __restrict  pmui,
                                          float * __restrict  Hzr,
                                          float * __restrict  Hzi ) {
                         
                         register __m512 H0r   = _mm512_loadu_ps(&pH0r[0]);
                         register __m512 H0i   = _mm512_loadu_ps(&pH0i[0]);
                         register __m512 k0z   = _mm512_loadu_ps(&pk0z[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 e0u0 = _mm512_set1_ps(0.00001763712109284471382861586f);
                         const __m512 pi4  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512 k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512 spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512 fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_ps(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_ps(k0z,sinps,_mm512_fmsub_ps(k0r,cosps,pi4));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_ps(k0r);
                         cmul_zmm16r4(H0r,H0i,cer,cei,&fracr,&fraci);
                         fracr = _mm512_mul_ps(fracr,cosps);
                         spsph = _mm512_mul_ps(sinps,sinph);
                         fraci = _mm512_mul_ps(fraci,cosps);
                         murp1 = _mm512_add_ps(mur,_1);
                         epsrp1= _mm512_add_ps(epsr,_1);
                         muip1 = _mm512_add_ps(mui,_1);
                         epsip1= _mm512_add_ps(epsi,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         t0r   = _mm512_mul_ps(e0u0,_mm512_div_ps(fracr,scpk0r));
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         t0i   = _mm512_mul_ps(e0u0,_mm512_div_ps(fraci,scpk0r));
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,spsph);
                         t0r  = _mm512_mul_ps(Ii,_mm512_mul_ps(t0r,k0a02));
                         divi = _mm512_mul_ps(divi,spsph);
                         t0i  = _mm512_mul_ps(Ii,_mm512_mul_ps(t0i,k0a02));
                         cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_storeu_ps(&Hzr[0], resr);
                         _mm512_storeu_ps(&Hzi[0], resi);
                 }


                  /*
                     Bistatic scattering width for (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                     Infinitely long cylinder.
                     TM-incident.
                     Formula 4.2-56
                 */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4256_zmm16r4(const __m512 a0,
                                            const __m512 k0a0,
                                            const __m512 psi,
                                            const __m512 phi,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {

                         const __m512 _4 = _mm512_set1_ps(4.0f);
                         const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512 epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512 epsrcps,epsicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512 t1r,t1i,cabs;
                         k0a03  = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrm1 = _mm512_sub_ps(epsr,_1);
                         spia   = _mm512_mul_ps(spi,a0);
                         epsim1 = _mm512_sub_ps(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_ps(epsr,_1);
                         cos2ps = _mm512_mul_ps(cosps,cosps);
                         epsip1 = _mm512_add_ps(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_ps(mur,_1);
                         sin2ps = _mm512_mul_ps(sinps,sinps);
                         muim1  = _mm512_sub_ps(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_ps(mur,_1);
                         t0     = _mm512_mul_ps(_4,cos2ps);
                         muip1  = _mm512_add_ps(mui,_1);
                         frac   = _mm512_div_ps(spia,t0);
                         epsrcps= _mm512_mul_ps(epsrm1,cos2ps);
                         frac   = _mm512_mul_ps(frac,k0a03);
                         epsicps= _mm512_mul_ps(epsim1,cos2ps);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_ps(_2,_mm512_mul_ps(divr,cosp));
                         t1r  = _mm512_sub_ps(epsrcps,t0r);
                         t0i  = _mm512_mul_ps(_2,_mm512_mul_ps(divi,cosp));
                         t1i  = _mm512_sub_ps(epsicps,t0i);
                         cabs = cabs_zmm16r4(t1r,t1i);
                         rcs  = _mm512_mul_ps(cabs,frac);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4256_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa0,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {

                         register __m512 a0   = _mm512_load_ps(&pa0[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 _4 = _mm512_set1_ps(4.0f);
                         const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512 epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512 epsrcps,epsicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512 t1r,t1i,cabs;
                         k0a03  = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrm1 = _mm512_sub_ps(epsr,_1);
                         spia   = _mm512_mul_ps(spi,a0);
                         epsim1 = _mm512_sub_ps(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_ps(epsr,_1);
                         cos2ps = _mm512_mul_ps(cosps,cosps);
                         epsip1 = _mm512_add_ps(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_ps(mur,_1);
                         sin2ps = _mm512_mul_ps(sinps,sinps);
                         muim1  = _mm512_sub_ps(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_ps(mur,_1);
                         t0     = _mm512_mul_ps(_4,cos2ps);
                         muip1  = _mm512_add_ps(mui,_1);
                         frac   = _mm512_div_ps(spia,t0);
                         epsrcps= _mm512_mul_ps(epsrm1,cos2ps);
                         frac   = _mm512_mul_ps(frac,k0a03);
                         epsicps= _mm512_mul_ps(epsim1,cos2ps);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_ps(_2,_mm512_mul_ps(divr,cosp));
                         t1r  = _mm512_sub_ps(epsrcps,t0r);
                         t0i  = _mm512_mul_ps(_2,_mm512_mul_ps(divi,cosp));
                         t1i  = _mm512_sub_ps(epsicps,t0i);
                         cabs = cabs_zmm16r4(t1r,t1i);
                         rcs  = _mm512_mul_ps(cabs,frac);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4256_zmm16r4_u(const float * __restrict  pa0,
                                              const float * __restrict  pk0a0,
                                              const float * __restrict  ppsi,
                                              const float * __restrict  pphi,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {

                         register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 _4 = _mm512_set1_ps(4.0f);
                         const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512 epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512 epsrcps,epsicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512 t1r,t1i,cabs;
                         k0a03  = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrm1 = _mm512_sub_ps(epsr,_1);
                         spia   = _mm512_mul_ps(spi,a0);
                         epsim1 = _mm512_sub_ps(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_ps(epsr,_1);
                         cos2ps = _mm512_mul_ps(cosps,cosps);
                         epsip1 = _mm512_add_ps(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_ps(mur,_1);
                         sin2ps = _mm512_mul_ps(sinps,sinps);
                         muim1  = _mm512_sub_ps(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_ps(mur,_1);
                         t0     = _mm512_mul_ps(_4,cos2ps);
                         muip1  = _mm512_add_ps(mui,_1);
                         frac   = _mm512_div_ps(spia,t0);
                         epsrcps= _mm512_mul_ps(epsrm1,cos2ps);
                         frac   = _mm512_mul_ps(frac,k0a03);
                         epsicps= _mm512_mul_ps(epsim1,cos2ps);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_ps(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_ps(mul2i,sin2ps,mul1i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_ps(_2,_mm512_mul_ps(divr,cosp));
                         t1r  = _mm512_sub_ps(epsrcps,t0r);
                         t0i  = _mm512_mul_ps(_2,_mm512_mul_ps(divi,cosp));
                         t1i  = _mm512_sub_ps(epsicps,t0i);
                         cabs = cabs_zmm16r4(t1r,t1i);
                         rcs  = _mm512_mul_ps(cabs,frac);
                         return (rcs);
                 }


                  /*
                     Bistatic scattering width for (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                     Infinitely long cylinder.
                     TE-incident.
                     Formula 4.2-58
                 */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4258_zmm16r4(const __m512 a0,
                                            const __m512 k0a0,
                                            const __m512 psi,
                                            const __m512 phi,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {

                         const __m512 _4 = _mm512_set1_ps(4.0f);
                         const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512 epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512 murcps,muicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512 t1r,t1i,cabs;
                         k0a03  = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrm1 = _mm512_sub_ps(epsr,_1);
                         spia   = _mm512_mul_ps(spi,a0);
                         epsim1 = _mm512_sub_ps(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_ps(epsr,_1);
                         cos2ps = _mm512_mul_ps(cosps,cosps);
                         epsip1 = _mm512_add_ps(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_ps(mur,_1);
                         sin2ps = _mm512_mul_ps(sinps,sinps);
                         muim1  = _mm512_sub_ps(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_ps(mur,_1);
                         t0     = _mm512_mul_ps(_4,cos2ps);
                         muip1  = _mm512_add_ps(mui,_1);
                         frac   = _mm512_div_ps(spia,t0);
                         murcps= _mm512_mul_ps(murm1,cos2ps);
                         frac   = _mm512_mul_ps(frac,k0a03);
                         muicps= _mm512_mul_ps(muim1,cos2ps);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_ps(mul1r,sin2ps,mul2r);
                         numi = _mm512_fmadd_ps(mul1i,sin2ps,mul2i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_ps(_2,_mm512_mul_ps(divr,cosp));
                         t1r  = _mm512_sub_ps(murcps,t0r);
                         t0i  = _mm512_mul_ps(_2,_mm512_mul_ps(divi,cosp));
                         t1i  = _mm512_sub_ps(muicps,t0i);
                         cabs = cabs_zmm16r4(t1r,t1i);
                         rcs  = _mm512_mul_ps(cabs,frac);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4258_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa0,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {

                         register __m512 a0    = _mm512_load_ps(&pa0[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 _4 = _mm512_set1_ps(4.0f);
                         const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512 epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512 murcps,muicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512 t1r,t1i,cabs;
                         k0a03  = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrm1 = _mm512_sub_ps(epsr,_1);
                         spia   = _mm512_mul_ps(spi,a0);
                         epsim1 = _mm512_sub_ps(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_ps(epsr,_1);
                         cos2ps = _mm512_mul_ps(cosps,cosps);
                         epsip1 = _mm512_add_ps(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_ps(mur,_1);
                         sin2ps = _mm512_mul_ps(sinps,sinps);
                         muim1  = _mm512_sub_ps(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_ps(mur,_1);
                         t0     = _mm512_mul_ps(_4,cos2ps);
                         muip1  = _mm512_add_ps(mui,_1);
                         frac   = _mm512_div_ps(spia,t0);
                         murcps= _mm512_mul_ps(murm1,cos2ps);
                         frac   = _mm512_mul_ps(frac,k0a03);
                         muicps= _mm512_mul_ps(muim1,cos2ps);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_ps(mul1r,sin2ps,mul2r);
                         numi = _mm512_fmadd_ps(mul1i,sin2ps,mul2i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_ps(_2,_mm512_mul_ps(divr,cosp));
                         t1r  = _mm512_sub_ps(murcps,t0r);
                         t0i  = _mm512_mul_ps(_2,_mm512_mul_ps(divi,cosp));
                         t1i  = _mm512_sub_ps(muicps,t0i);
                         cabs = cabs_zmm16r4(t1r,t1i);
                         rcs  = _mm512_mul_ps(cabs,frac);
                         return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4258_zmm16r4_u(const float * __restrict  pa0,
                                              const float * __restrict  pk0a0,
                                              const float * __restrict  ppsi,
                                              const float * __restrict  pphi,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {

                         register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 _4 = _mm512_set1_ps(4.0f);
                         const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2 = _mm512_set1_ps(2.0f);
                         register __m512 rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512 epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512 murcps,muicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512 t1r,t1i,cabs;
                         k0a03  = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrm1 = _mm512_sub_ps(epsr,_1);
                         spia   = _mm512_mul_ps(spi,a0);
                         epsim1 = _mm512_sub_ps(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_ps(epsr,_1);
                         cos2ps = _mm512_mul_ps(cosps,cosps);
                         epsip1 = _mm512_add_ps(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_ps(mur,_1);
                         sin2ps = _mm512_mul_ps(sinps,sinps);
                         muim1  = _mm512_sub_ps(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_ps(mur,_1);
                         t0     = _mm512_mul_ps(_4,cos2ps);
                         muip1  = _mm512_add_ps(mui,_1);
                         frac   = _mm512_div_ps(spia,t0);
                         murcps= _mm512_mul_ps(murm1,cos2ps);
                         frac   = _mm512_mul_ps(frac,k0a03);
                         muicps= _mm512_mul_ps(muim1,cos2ps);
                         cmul_zmm16r4(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm16r4(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_ps(mul1r,sin2ps,mul2r);
                         numi = _mm512_fmadd_ps(mul1i,sin2ps,mul2i);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm16r4(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_ps(_2,_mm512_mul_ps(divr,cosp));
                         t1r  = _mm512_sub_ps(murcps,t0r);
                         t0i  = _mm512_mul_ps(_2,_mm512_mul_ps(divi,cosp));
                         t1i  = _mm512_sub_ps(muicps,t0i);
                         cabs = cabs_zmm16r4(t1r,t1i);
                         rcs  = _mm512_mul_ps(cabs,frac);
                         return (rcs);
                 }


                   /*

                           Bistatic scattering width for (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                           Infinitely long cylinder.
                           TM-incident.
                           Formula 4.2-57
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4257_zmm16r4(const __m512 a0,
                                            const __m512 k0a0,
                                            const __m512 psi,
                                            const __m512 phi,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {

                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 spi = _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2  = _mm512_set1_ps(2.0f);
                         register __m512 spi4,rcs,cos2ps,sinps,sinp,k0a03;
                         register __m512 frac,divr,divi,mul1r,mul1i,mul2r,mul2i;
                         register __m512 epsrp1,epsip1,murp1,muip1,t0,cabs;
                         k0a03 = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrp1= _mm512_add_ps(epsr,_1);
                         cosps = xcosf(psi);
                         epsip1= _mm512_add_ps(epsi,_1);
                         spi4  = _mm512_mul_ps(spi,_mm512_mul_ps(a0,a0));
                         spi4  = _mm512_mul_ps(_4,spi4);
                         murp1 = _mm512_add_ps(mur,_1);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         muip1 = _mm512_add_ps(mui,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         sinps = xsinf(psi);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         sinp  = xsinf(phi);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         frac  = _mm512_mul_ps(k0a03,_mm512_div_ps(spi4,cos2ps));
                         t0    = _mm512_mul_ps(sinps,sinp);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,t0);
                         divi = _mm512_mul_ps(divi,t0);
                         cabs = cabs_zmm16r4(divr,divi);
                         rcs  = _mm512_mul_ps(frac,cabs);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4257_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa0,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a0,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {

                         register __m512 a0    = _mm512_load_ps(&pa0[0]);
                         register __m512 k0a0  = _mm512_load_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_load_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_load_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_load_ps(&pmur[0]);
                         register __m512 mui   = _mm512_load_ps(&pmui[0]); 
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 spi = _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2  = _mm512_set1_ps(2.0f);
                         register __m512 spi4,rcs,cos2ps,sinps,sinp,k0a03;
                         register __m512 frac,divr,divi,mul1r,mul1i,mul2r,mul2i;
                         register __m512 epsrp1,epsip1,murp1,muip1,t0,cabs;
                         k0a03 = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrp1= _mm512_add_ps(epsr,_1);
                         cosps = xcosf(psi);
                         epsip1= _mm512_add_ps(epsi,_1);
                         spi4  = _mm512_mul_ps(spi,_mm512_mul_ps(a0,a0));
                         spi4  = _mm512_mul_ps(_4,spi4);
                         murp1 = _mm512_add_ps(mur,_1);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         muip1 = _mm512_add_ps(mui,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         sinps = xsinf(psi);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         sinp  = xsinf(phi);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         frac  = _mm512_mul_ps(k0a03,_mm512_div_ps(spi4,cos2ps));
                         t0    = _mm512_mul_ps(sinps,sinp);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,t0);
                         divi = _mm512_mul_ps(divi,t0);
                         cabs = cabs_zmm16r4(divr,divi);
                         rcs  = _mm512_mul_ps(frac,cabs);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4257_zmm16r4_u(const float * __restrict  pa0,
                                              const float * __restrict  pk0a0,
                                              const float * __restrict  ppsi,
                                              const float * __restrict  pphi,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {

                         register __m512 a0    = _mm512_loadu_ps(&pa0[0]);
                         register __m512 k0a0  = _mm512_loadu_ps(&pk0a0[0]);
                         register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 pphi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                         register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                         register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                         register __m512 mui   = _mm512_loadu_ps(&pmui[0]); 
                         const __m512 _4  = _mm512_set1_ps(4.0f);
                         const __m512 spi = _mm512_set1_ps(9.869604401089358618834490999876f);
                         const __m512 _2  = _mm512_set1_ps(2.0f);
                         register __m512 spi4,rcs,cos2ps,sinps,sinp,k0a03;
                         register __m512 frac,divr,divi,mul1r,mul1i,mul2r,mul2i;
                         register __m512 epsrp1,epsip1,murp1,muip1,t0,cabs;
                         k0a03 = _mm512_mul_ps(k0a0,_mm512_mul_ps(k0a0,k0a0));
                         epsrp1= _mm512_add_ps(epsr,_1);
                         cosps = xcosf(psi);
                         epsip1= _mm512_add_ps(epsi,_1);
                         spi4  = _mm512_mul_ps(spi,_mm512_mul_ps(a0,a0));
                         spi4  = _mm512_mul_ps(_4,spi4);
                         murp1 = _mm512_add_ps(mur,_1);
                         cos2ps= _mm512_mul_ps(cosps,cosps);
                         muip1 = _mm512_add_ps(mui,_1);
                         cmul_zmm16r4(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         sinps = xsinf(psi);
                         mul1r = _mm512_sub_ps(mul1r,_1);
                         mul1i = _mm512_sub_ps(mul1i,_1);
                         sinp  = xsinf(phi);
                         cmul_zmm16r4(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         frac  = _mm512_mul_ps(k0a03,_mm512_div_ps(spi4,cos2ps));
                         t0    = _mm512_mul_ps(sinps,sinp);
                         cdiv_zmm16r4(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_ps(divr,t0);
                         divi = _mm512_mul_ps(divi,t0);
                         cabs = cabs_zmm16r4(divr,divi);
                         rcs  = _mm512_mul_ps(frac,cabs);
                         return (rcs);
                 }


                   /*
                       Circular cylinders of finite length.
                       Cylinder radius small (k0a<1.0)
                       Wire limit of cylinder (h>>a).
                       E-field
                       Formula 4.3-9
                   */

                    
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f439_zmm16r4(const __m512 EIr,
                                        const __m512 EIi,
                                        const __m512 r,
                                        const __m512 k0,
                                        const __m512 psii,
                                        const __m512 psis,
                                        const __m512 h,
                                        const __m512 ln4ha,
                                        __m512 * __restrict ESr,
                                        __m512 * __restrict ESi) {

                       const __m512 thrd = _mm512_set1_ps(0.333333333333333333333333333333f);
                       const __m512 _1   = _mm512_set1_ps(1.0f);
                       register __m512 ir,k02,ear,eai,cer,cei,h3,cpsii,cpsis;
                       register __m512 num,rat,den,t0r,t0i,mulr,muli;
                       cpsii = xcosf(psii);
                       k02   = _mm512_mul_ps(thrd,_mm512_mul_ps(k0,k0));
                       ir    = _mm512_rcp14_ps(r);
                       ear   = Ir;
                       cpsis = xcosf(psis);
                       eai   = _mm512_mul_ps(k0,r);
                       den   = _mm512_sub_ps(ln4ha,_1);
                       h3    = _mm512_mul_ps(h,_mm512_mul_ps(h,h));
                       cexp_zmm16r4(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_ps(cer,ir);
                       num   = _mm512_mul_ps(h3,_mm512_mul_ps(cpsis,cpsii));
                       cei   = _mm512_mul_ps(cei,ir);
                       rat   = _mm512_div_ps(num,den);
                       t0r   = _mm512_mul_ps(EIr,rat);
                       t0i   = _mm512_mul_ps(EIi,rat);
                       cmul_zmm16r4(cer,cei,t0r,t0i,&mulr,&muli);
                       *ESr = _mm512_mul_ps(mulr,k02);
                       *ESi = _mm512_mul_ps(muli,k02);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f439_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEIr,
                                          const float * __restrict __ATTR_ALIGN__(64) pEIi,
                                          const float * __restrict __ATTR_ALIGN__(64) pr,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                          const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                          const float * __restrict __ATTR_ALIGN__(64) ph,
                                          const float * __restrict __ATTR_ALIGN__(64) pln4ha,
                                          float * __restrict __ATTR_ALIGN__(64) ESr,
                                          float * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512 EIr  = _mm512_load_ps(&pEIr[0]);
                       register __m512 EIi  = _mm512_load_ps(&pEIi[0]);
                       register __m512 r    = _mm512_load_ps(&pr[0]);
                       register __m512 k0   = _mm512_load_ps(&pk0[0]);
                       register __m512 psii = _mm512_load_ps(&ppsii[0]);
                       register __m512 psis = _mm512_load_ps(&ppsis[0]);
                       register __m512 h    = _mm512_load_ps(&ph[0]);
                       register __m512 ln4ha= _mm512_load_ps(&pln4ha[0]);
                       const __m512 thrd = _mm512_set1_ps(0.333333333333333333333333333333f);
                       const __m512 _1   = _mm512_set1_ps(1.0f);
                       register __m512 ir,k02,ear,eai,cer,cei,h3,cpsii,cpsis;
                       register __m512 num,rat,den,t0r,t0i,mulr,muli;
                       cpsii = xcosf(psii);
                       k02   = _mm512_mul_ps(thrd,_mm512_mul_ps(k0,k0));
                       ir    = _mm512_rcp14_ps(r);
                       ear   = Ir;
                       cpsis = xcosf(psis);
                       eai   = _mm512_mul_ps(k0,r);
                       den   = _mm512_sub_ps(ln4ha,_1);
                       h3    = _mm512_mul_ps(h,_mm512_mul_ps(h,h));
                       cexp_zmm16r4(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_ps(cer,ir);
                       num   = _mm512_mul_ps(h3,_mm512_mul_ps(cpsis,cpsii));
                       cei   = _mm512_mul_ps(cei,ir);
                       rat   = _mm512_div_ps(num,den);
                       t0r   = _mm512_mul_ps(EIr,rat);
                       t0i   = _mm512_mul_ps(EIi,rat);
                       cmul_zmm16r4(cer,cei,t0r,t0i,&mulr,&muli);
                       _mm512_store_ps(&ESr[0] ,_mm512_mul_ps(mulr,k02));
                       _mm512_store_ps(&ESi[0] ,_mm512_mul_ps(muli,k02));
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f439_zmm16r4_u(const float * __restrict  pEIr,
                                          const float * __restrict  pEIi,
                                          const float * __restrict  pr,
                                          const float * __restrict  pk0,
                                          const float * __restrict  ppsii,
                                          const float * __restrict  ppsis,
                                          const float * __restrict  ph,
                                          const float * __restrict  pln4ha,
                                          float * __restrict  ESr,
                                          float * __restrict  ESi) {

                       register __m512 EIr  = _mm512_loadu_ps(&pEIr[0]);
                       register __m512 EIi  = _mm512_loadu_ps(&pEIi[0]);
                       register __m512 r    = _mm512_loadu_ps(&pr[0]);
                       register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                       register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                       register __m512 psis = _mm512_loadu_ps(&ppsis[0]);
                       register __m512 h    = _mm512_loadu_ps(&ph[0]);
                       register __m512 ln4ha= _mm512_loadu_ps(&pln4ha[0]);
                       const __m512 thrd = _mm512_set1_ps(0.333333333333333333333333333333f);
                       const __m512 _1   = _mm512_set1_ps(1.0f);
                       register __m512 ir,k02,ear,eai,cer,cei,h3,cpsii,cpsis;
                       register __m512 num,rat,den,t0r,t0i,mulr,muli;
                       cpsii = xcosf(psii);
                       k02   = _mm512_mul_ps(thrd,_mm512_mul_ps(k0,k0));
                       ir    = _mm512_rcp14_ps(r);
                       ear   = Ir;
                       cpsis = xcosf(psis);
                       eai   = _mm512_mul_ps(k0,r);
                       den   = _mm512_sub_ps(ln4ha,_1);
                       h3    = _mm512_mul_ps(h,_mm512_mul_ps(h,h));
                       cexp_zmm16r4(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_ps(cer,ir);
                       num   = _mm512_mul_ps(h3,_mm512_mul_ps(cpsis,cpsii));
                       cei   = _mm512_mul_ps(cei,ir);
                       rat   = _mm512_div_ps(num,den);
                       t0r   = _mm512_mul_ps(EIr,rat);
                       t0i   = _mm512_mul_ps(EIi,rat);
                       cmul_zmm16r4(cer,cei,t0r,t0i,&mulr,&muli);
                       _mm512_storeu_ps(&ESr[0] ,_mm512_mul_ps(mulr,k02));
                       _mm512_storeu_ps(&ESi[0] ,_mm512_mul_ps(muli,k02));
                 }


                  /*
                       Circular cylinders of finite length.
                       Cylinder radius small (k0a<1.0)
                       Wire limit of cylinder (h>>a).
                       RCS.
                       Formula 4.3-10

                    */
                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4310_zmm16r4(const __m512 k0,
                                            const __m512 h,
                                            const __m512 psii,
                                            const __m512 psis,
                                            const __m512 ln4h) {

                          const __m512 _4pi9 = _mm512_set1_ps(1.396263401595463661538952614791f);
                          const __m512 _1    = _mm512_set1_ps(1.0f);
                          register __m512 cpsii,cpsis,c2psii,c2psis,den,num,t0;
                          register __m512 k04,h6,rcs,t1,h2,rat,frac;
                          h2     = _mm512_mul_ps(h,h);
                          k04    = _mm512_mul_ps(_mm512_mul_ps(k0,k0),
                                              _mm512_mul_ps(k0,k0));
                          t0     = _mm512_sub_ps(ln4h,_1);
                          cpsii  = xcosf(psii);
                          c2psii = _mm512_mul_ps(cpsii,cpsii);
                          den    = _mm512_mul_ps(t0,t0); 
                          t1     = _mm512_mul_ps(h,h2);
                          h6     = _mm512_mul_ps(t1,h2);
                          cpsis  = xcosf(psis);
                          c2psis = _mm512_mul_ps(cpsis,cpsis);
                          num    = _mm512_mul_ps(c2psis,c2psii);
                          frac   = _mm512_mul_ps(_4pi9,_mm512_mul_ps(k04,h6));
                          rat    = _mm512_div_ps(num,den);
                          rcs    = _mm512_mul_ps(frac,rat);
                          return (rcs);
                 } 


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4310_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const float * __restrict __ATTR_ALIGN__(64)  ph,
                                              const float * __restrict __ATTR_ALIGN__(64)  ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64)  ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64)  pln4h) {

                          register __m512 k0   = _mm512_load_ps(&pk0[0]);
                          register __m512 h    = _mm512_load_ps(&ph[0]);
                          register __m512 psii = _mm512_load_ps(&ppsii[0]);
                          register __m512 psis = _mm512_load_ps(&ppsis[0]);
                          register __m512 pln4h= _mm512_load_ps(&pln4h[0]);
                          const __m512 _4pi9 = _mm512_set1_ps(1.396263401595463661538952614791f);
                          const __m512 _1    = _mm512_set1_ps(1.0f);
                          register __m512 cpsii,cpsis,c2psii,c2psis,den,num,t0;
                          register __m512 k04,h6,rcs,t1,h2,rat,frac;
                          h2     = _mm512_mul_ps(h,h);
                          k04    = _mm512_mul_ps(_mm512_mul_ps(k0,k0),
                                              _mm512_mul_ps(k0,k0));
                          t0     = _mm512_sub_ps(ln4h,_1);
                          cpsii  = xcosf(psii);
                          c2psii = _mm512_mul_ps(cpsii,cpsii);
                          den    = _mm512_mul_ps(t0,t0); 
                          t1     = _mm512_mul_ps(h,h2);
                          h6     = _mm512_mul_ps(t1,h2);
                          cpsis  = xcosf(psis);
                          c2psis = _mm512_mul_ps(cpsis,cpsis);
                          num    = _mm512_mul_ps(c2psis,c2psii);
                          frac   = _mm512_mul_ps(_4pi9,_mm512_mul_ps(k04,h6));
                          rat    = _mm512_div_ps(num,den);
                          rcs    = _mm512_mul_ps(frac,rat);
                          return (rcs);
                 } 


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4310_zmm16r4_u(const float * __restrict   pk0,
                                              const float * __restrict   ph,
                                              const float * __restrict   ppsii,
                                              const float * __restrict   ppsis,
                                              const float * __restrict   pln4h) {

                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                          register __m512 h    = _mm512_loadu_ps(&ph[0]);
                          register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                          register __m512 psis = _mm512_loadu_ps(&ppsis[0]);
                          register __m512 pln4h= _mm512_loadu_ps(&pln4h[0]);
                          const __m512 _4pi9 = _mm512_set1_ps(1.396263401595463661538952614791f);
                          const __m512 _1    = _mm512_set1_ps(1.0f);
                          register __m512 cpsii,cpsis,c2psii,c2psis,den,num,t0;
                          register __m512 k04,h6,rcs,t1,h2,rat,frac;
                          h2     = _mm512_mul_ps(h,h);
                          k04    = _mm512_mul_ps(_mm512_mul_ps(k0,k0),
                                              _mm512_mul_ps(k0,k0));
                          t0     = _mm512_sub_ps(ln4h,_1);
                          cpsii  = xcosf(psii);
                          c2psii = _mm512_mul_ps(cpsii,cpsii);
                          den    = _mm512_mul_ps(t0,t0); 
                          t1     = _mm512_mul_ps(h,h2);
                          h6     = _mm512_mul_ps(t1,h2);
                          cpsis  = xcosf(psis);
                          c2psis = _mm512_mul_ps(cpsis,cpsis);
                          num    = _mm512_mul_ps(c2psis,c2psii);
                          frac   = _mm512_mul_ps(_4pi9,_mm512_mul_ps(k04,h6));
                          rat    = _mm512_div_ps(num,den);
                          rcs    = _mm512_mul_ps(frac,rat);
                          return (rcs);
                 } 


                  /*
                         The average dipole scattering RCS when the incidence
                         and scattered polarization direction coincide.
                         Formula 4.3-11
                    */

                   
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4311_zmm16r4(const __m512 k0,
                                            const __m512 h,
                                            const __m512 ln4h) {

                          const __m512 _4pi45 = _mm512_set1_ps(0.279252680319092732307790522958f);
                          const __m512 _1     = _mm512_set1_ps(1.0f);
                          register __m512 rcs,den,inv,k04,h6,h2,t0,t1;
                          h2     = _mm512_mul_ps(h,h);
                          k04    = _mm512_mul_ps(_mm512_mul_ps(k0,k0),
                                              _mm512_mul_ps(k0,k0));
                          t0     = _mm512_sub_ps(ln4h,_1);
                          t1     = _mm512_mul_ps(h,h2);
                          h6     = _mm512_mul_ps(t1,h2);
                          den    = _mm512_mul_ps(t0,t0);
                          inv    = _mm512_div_ps(_1,den);
                          t0     = _mm512_mul_ps(_4pi45,_mm512_mul_ps(k04,h6));
                          rcs    = _mm512_mul_ps(t0,inv);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4311_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) ph,
                                              const float * __restrict __ATTR_ALIGN__(64) pln4h) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 h   = _mm512_load_ps(&ph[0]);
                          register __m512 ln4h= _mm512_load_ps(&pln4h[0]);
                          const __m512 _4pi45 = _mm512_set1_ps(0.279252680319092732307790522958f);
                          const __m512 _1     = _mm512_set1_ps(1.0f);
                          register __m512 rcs,den,inv,k04,h6,h2,t0,t1;
                          h2     = _mm512_mul_ps(h,h);
                          k04    = _mm512_mul_ps(_mm512_mul_ps(k0,k0),
                                              _mm512_mul_ps(k0,k0));
                          t0     = _mm512_sub_ps(ln4h,_1);
                          t1     = _mm512_mul_ps(h,h2);
                          h6     = _mm512_mul_ps(t1,h2);
                          den    = _mm512_mul_ps(t0,t0);
                          inv    = _mm512_div_ps(_1,den);
                          t0     = _mm512_mul_ps(_4pi45,_mm512_mul_ps(k04,h6));
                          rcs    = _mm512_mul_ps(t0,inv);
                          return (rcs);
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4311_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  ph,
                                              const float * __restrict  pln4h) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 h   = _mm512_loadu_ps(&ph[0]);
                          register __m512 ln4h= _mm512_loadu_ps(&pln4h[0]);
                          const __m512 _4pi45 = _mm512_set1_ps(0.279252680319092732307790522958f);
                          const __m512 _1     = _mm512_set1_ps(1.0f);
                          register __m512 rcs,den,inv,k04,h6,h2,t0,t1;
                          h2     = _mm512_mul_ps(h,h);
                          k04    = _mm512_mul_ps(_mm512_mul_ps(k0,k0),
                                              _mm512_mul_ps(k0,k0));
                          t0     = _mm512_sub_ps(ln4h,_1);
                          t1     = _mm512_mul_ps(h,h2);
                          h6     = _mm512_mul_ps(t1,h2);
                          den    = _mm512_mul_ps(t0,t0);
                          inv    = _mm512_div_ps(_1,den);
                          t0     = _mm512_mul_ps(_4pi45,_mm512_mul_ps(k04,h6));
                          rcs    = _mm512_mul_ps(t0,inv);
                          return (rcs);
               }


                  /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-18
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4318_zmm16r4(const __m512 EIr,
                                         const __m512 EIi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 psii,
                                         const __m512 psis,
                                         const __m512 phi,
                                         const __m512 a,
                                         __m512 * __restrict ESr,
                                         __m512 * __restrict ESi) {

                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         register __m512 ir,a3,k02,cosp,spsii,spsis,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         spsis= xsinf(psis);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_mul_ps(spsii,_mm512_mul_ps(spsis,cosp));
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4318_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  ppsii,
                                           const float * __restrict __ATTR_ALIGN__(64)  ppsis,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pa,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64)  ESi) {

                         register __m512 EIr  = _mm512_load_ps(&pEIr[0]);
                         register __m512 EIi  = _mm512_load_ps(&pEIi[0]);
                         register __m512 k0   = _mm512_load_ps(&pk0[0]);
                         register __m512 r    = _mm512_load_ps(&pr[0]);
                         register __m512 psii = _mm512_load_ps(&ppsii[0]);
                         register __m512 psis = _mm512_load_ps(&ppsis[0]);
                         register __m512 phi  = _mm512_load_ps(&pphi[0]);
                         register __m512 a    = _mm512_load_ps(&a[0]);
                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         register __m512 ir,a3,k02,cosp,spsii,spsis,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         spsis= xsinf(psis);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_mul_ps(spsii,_mm512_mul_ps(spsis,cosp));
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_store_ps(&ESr[0], resr);
                         _mm512_store_ps(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4318_zmm16r4_u(const float * __restrict   pEIr,
                                           const float * __restrict   pEIi,
                                           const float * __restrict   pk0,
                                           const float * __restrict   pr,
                                           const float * __restrict   ppsii,
                                           const float * __restrict   ppsis,
                                           const float * __restrict   pphi,
                                           const float * __restrict   pa,
                                           float * __restrict  ESr,
                                           float * __restrict   ESi) {

                         register __m512 EIr  = _mm512_loadu_ps(&pEIr[0]);
                         register __m512 EIi  = _mm512_loadu_ps(&pEIi[0]);
                         register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                         register __m512 r    = _mm512_loadu_ps(&pr[0]);
                         register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                         register __m512 psis = _mm512_loadu_ps(&ppsis[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 a    = _mm512_loadu_ps(&a[0]);
                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         register __m512 ir,a3,k02,cosp,spsii,spsis,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         spsis= xsinf(psis);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_mul_ps(spsii,_mm512_mul_ps(spsis,cosp));
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_storeu_ps(&ESr[0], resr);
                         _mm512_storeu_ps(&ESi[0], resi);
                }


                   /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-19
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4319_zmm16r4(const __m512 EIr,
                                         const __m512 EIi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 psii,
                                         const __m512 phi,
                                         const __m512 a,
                                         __m512 * __restrict ESr,
                                         __m512 * __restrict ESi) {

                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         register __m512 ir,a3,k02,sinp,spsii,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         sinp = xsinf(phi);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_mul_ps(spsii,sinp);
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4319_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  ppsii,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pa,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64)  ESi) {

                         register __m512 EIr  = _mm512_load_ps(&pEIr[0]);
                         register __m512 EIi  = _mm512_load_ps(&pEIi[0]);
                         register __m512 k0   = _mm512_load_ps(&pk0[0]);
                         register __m512 r    = _mm512_load_ps(&pr[0]);
                         register __m512 psii = _mm512_load_ps(&ppsii[0]);
                         register __m512 phi  = _mm512_load_ps(&pphi[0]);
                         register __m512 a    = _mm512_load_ps(&a[0]);
                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         register __m512 ir,a3,k02,sinp,spsii,t0,t1,resr,resi;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         sinp = xsinf(phi);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_mul_ps(spsii,sinp);
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_store_ps(&ESr[0], resr);
                         _mm512_store_ps(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4319_zmm16r4_u(const float * __restrict   pEIr,
                                           const float * __restrict   pEIi,
                                           const float * __restrict   pk0,
                                           const float * __restrict   pr,
                                           const float * __restrict   ppsii,
                                           const float * __restrict   pphi,
                                           const float * __restrict   pa,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                         register __m512 EIr  = _mm512_loadu_ps(&pEIr[0]);
                         register __m512 EIi  = _mm512_loadu_ps(&pEIi[0]);
                         register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                         register __m512 r    = _mm512_loadu_ps(&pr[0]);
                         register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 a    = _mm512_loadu_ps(&a[0]);
                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         register __m512 ir,a3,k02,sinp,spsii,t0,t1,resr,resi;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         sinp = xsinf(phi);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_mul_ps(spsii,sinp);
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_storeu_ps(&ESr[0], resr);
                         _mm512_storeu_ps(&ESi[0], resi);
                }


                 /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-20
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4320_zmm16r4(const __m512 EIr,
                                         const __m512 EIi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 psis,
                                         const __m512 phi,
                                         const __m512 a,
                                         __m512 * __restrict ESr,
                                         __m512 * __restrict ESi) {
                      
                        ES_f4319_zmm16r4(EIr,EIi,k0,r,psis,phi,a,&ESr,&ESi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4320_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  ppsis,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pa,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64)  ESi) {

                      ES_f4319_zmm16r4_a(pEIr,pEIi,pk0,pr,ppsis,pphi,ESr,ESi);
              }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4320_zmm16r4_u(const float * __restrict   pEIr,
                                           const float * __restrict   pEIi,
                                           const float * __restrict   pk0,
                                           const float * __restrict   pr,
                                           const float * __restrict   ppsis,
                                           const float * __restrict   pphi,
                                           const float * __restrict   pa,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                      ES_f4319_zmm16r4_u(pEIr,pEIi,pk0,pr,ppsis,pphi,ESr,ESi);
              }


                   /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-21
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4321_zmm16r4(const __m512 EIr,
                                         const __m512 EIi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 psii,
                                         const __m512 psis,
                                         const __m512 phi,
                                         const __m512 a,
                                         __m512 * __restrict ESr,
                                         __m512 * __restrict ESi) {

                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         const __m512 hlf   = _mm512_set1_ps(0.5f);
                         register __m512 ir,a3,k02,cosp,cpsii,cpsis,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         cpsis= xcosf(psis);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         cpsii= xcosf(psii);
                         cpsii= _mm512_mul_ps(hlf,cpsii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_fmadd_ps(cpsii,cpsis,cosp);
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,*ESr,*ESi);
              }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4321_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  ppsii,
                                           const float * __restrict __ATTR_ALIGN__(64)  ppsis,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pa,
                                           float * __restrict __ATTR_ALIGN__(64) ESr,
                                           float * __restrict __ATTR_ALIGN__(64)  ESi) {

                         register __m512 EIr  = _mm512_load_ps(&pEIr[0]);
                         register __m512 EIi  = _mm512_load_ps(&pEIi[0]);
                         register __m512 k0   = _mm512_load_ps(&pk0[0]);
                         register __m512 r    = _mm512_load_ps(&pr[0]);
                         register __m512 psii = _mm512_load_ps(&ppsii[0]);
                         register __m512 psis = _mm512_load_ps(&ppsis[0]);
                         register __m512 phi  = _mm512_load_ps(&pphi[0]);
                         register __m512 a    = _mm512_load_ps(&a[0]);
                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         const __m512 hlf   = _mm512_set1_ps(0.5f);
                         register __m512 ir,a3,k02,cosp,cpsii,cpsis,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         cpsis= xcosf(psis);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         cpsii= xcosf(psii);
                         cpsii= _mm512_mul_ps(hlf,cpsii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_fmadd_ps(cpsii,cpsis,cosp);
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_store_ps(&ESr[0], resr);
                         _mm512_store_ps(&ESi[0], resi);
              }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4321_zmm16r4_u(const float * __restrict   pEIr,
                                           const float * __restrict   pEIi,
                                           const float * __restrict   pk0,
                                           const float * __restrict   pr,
                                           const float * __restrict   ppsii,
                                           const float * __restrict   ppsis,
                                           const float * __restrict   pphi,
                                           const float * __restrict   pa,
                                           float * __restrict  ESr,
                                           float * __restrict  ESi) {

                         register __m512 EIr  = _mm512_loadu_ps(&pEIr[0]);
                         register __m512 EIi  = _mm512_loadu_ps(&pEIi[0]);
                         register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                         register __m512 r    = _mm512_loadu_ps(&pr[0]);
                         register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                         register __m512 psis = _mm512_loadu_ps(&ppsis[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                         register __m512 a    = _mm512_loadu_ps(&a[0]);
                         const __m512 _43pi = _mm512_set1_ps(0.424413181578387562050356702327f);
                         const __m512 hlf   = _mm512_set1_ps(0.5f);
                         register __m512 ir,a3,k02,cosp,cpsii,cpsis,t0,t1;
                         register __m512 ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                         ir   = _mm512_rcp14_ps(r);
                         cpsis= xcosf(psis);
                         k02  = _mm512_mul_ps(k0,k0);
                         ear  = _mm512_mul_ps(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_ps(_43pi,k02);
                         cpsii= xcosf(psii);
                         cpsii= _mm512_mul_ps(hlf,cpsii);
                         t0r  = _mm512_mul_ps(t0,_mm512_mul_ps(cer,ir));
                         t0i  = _mm512_mul_ps(t0,_mm512_mul_ps(cei,ir)); 
                         t1   = _mm512_fmadd_ps(cpsii,cpsis,cosp);
                         t0   = _mm512_mul_ps(a3,t1);
                         t1r  = _mm512_mul_ps(EIr,t0);
                         t1i  = _mm512_mul_ps(EIi,t0);
                         cmul_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_storeu_ps(&ESr[0], resr);
                         _mm512_storeu_ps(&ESi[0], resi);
              }


                 /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-22
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4322_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 psii,
                                            const __m512 psis,
                                            const __m512 phi) {

                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 rcs,k04,a6,t0,t1,spsii,spsis,cosp;
                          register __m512 s2psii,s2psis,cos2p,t2;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          t1    = _mm512_mul_ps(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          s2psii= _mm512_mul_ps(spsii,spsii);
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          spsis = xsinf(psis);
                          s2psis= _mm512_mul_ps(psis,psis);
                          t3    = _mm512_mul_ps(s2psii,_mm512_mul_ps(s2psis,cosp));
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4322_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512  k0  = _mm512_load_ps(&pk0[0]);
                          register __m512  a   = _mm512_load_ps(&pa[0]);
                          register __m512  psii= _mm512_load_ps(&ppsii[0]);
                          register __m512  psis= _mm512_load_ps(&ppsis[0]);
                          register __m512  phi = _mm512_load_ps(&pphi[0]);
                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 rcs,k04,a6,t0,t1,spsii,spsis,cosp;
                          register __m512 s2psii,s2psis,cos2p,t2;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          t1    = _mm512_mul_ps(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          s2psii= _mm512_mul_ps(spsii,spsii);
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          spsis = xsinf(psis);
                          s2psis= _mm512_mul_ps(psis,psis);
                          t3    = _mm512_mul_ps(s2psii,_mm512_mul_ps(s2psis,cosp));
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4322_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  ppsii,
                                              const float * __restrict  ppsis,
                                              const float * __restrict  pphi) {

                          register __m512  k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512  a   = _mm512_loadu_ps(&pa[0]);
                          register __m512  psii= _mm512_loadu_ps(&ppsii[0]);
                          register __m512  psis= _mm512_loadu_ps(&ppsis[0]);
                          register __m512  phi = _mm512_loadu_ps(&pphi[0]);
                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 rcs,k04,a6,t0,t1,spsii,spsis,cosp;
                          register __m512 s2psii,s2psis,cos2p,t2;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          cos2p = _mm512_mul_ps(cosp,cosp);
                          t1    = _mm512_mul_ps(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          s2psii= _mm512_mul_ps(spsii,spsii);
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          spsis = xsinf(psis);
                          s2psis= _mm512_mul_ps(psis,psis);
                          t3    = _mm512_mul_ps(s2psii,_mm512_mul_ps(s2psis,cosp));
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,t3));
                          return (rcs);
                }


                    /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-23
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4323_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 psii,
                                            const __m512 phi) {

                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 rcs,k04,a6,t0,t1,spsii,sinp;
                          register __m512 s2psii,sin2p,t2;
                          sinp  = xsinf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          sin2p = _mm512_mul_ps(sinp,sinp);
                          t1    = _mm512_mul_ps(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          s2psii= _mm512_mul_ps(spsii,spsii);
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          t3    = _mm512_mul_ps(s2psii,sin2p);
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4323_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512  k0  = _mm512_load_ps(&pk0[0]);
                          register __m512  a   = _mm512_load_ps(&pa[0]);
                          register __m512  psii= _mm512_load_ps(&ppsii[0]);
                          register __m512  phi = _mm512_load_ps(&pphi[0]);
                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 rcs,k04,a6,t0,t1,spsii,sinp;
                          register __m512 s2psii,sin2p,t2;
                          sinp  = xsinf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          sin2p = _mm512_mul_ps(sinp,sinp);
                          t1    = _mm512_mul_ps(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          s2psii= _mm512_mul_ps(spsii,spsii);
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          t3    = _mm512_mul_ps(s2psii,sin2p);
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4323_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  ppsii,
                                              const float * __restrict  pphi) {

                          register __m512  k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512  a   = _mm512_loadu_ps(&pa[0]);
                          register __m512  psii= _mm512_loadu_ps(&ppsii[0]);
                          register __m512  phi = _mm512_loadu_ps(&pphi[0]);
                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 rcs,k04,a6,t0,t1,spsii,sinp;
                          register __m512 s2psii,sin2p,t2;
                          sinp  = xsinf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          sin2p = _mm512_mul_ps(sinp,sinp);
                          t1    = _mm512_mul_ps(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          s2psii= _mm512_mul_ps(spsii,spsii);
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          t3    = _mm512_mul_ps(s2psii,sin2p);
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,t3));
                          return (rcs);
                }


                  /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-24
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4324_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 psis,
                                            const __m512 phi) {

                          register __m512 rcs;
                          rcs = rcs_f4323_zmm16r4(k0,a,psis,phi);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4324_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {

                           register __m512 rcs;
                           rcs = rcs_f4323_zmm16r4_a(pk0,pa,ppsis,pphi);
                           return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4324_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  ppsis,
                                              const float * __restrict  pphi) {

                           register __m512 rcs;
                           rcs = rcs_f4323_zmm16r4_u(pk0,pa,ppsis,pphi);
                           return (rcs);
                }


                  /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-25
                   */

  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4325_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 psii,
                                            const __m512 psis,
                                            const __m512 phi) {

                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          const __m512 hlf    = _mm512_set1_ps(0.5f);
                          register __m512 rcs,k04,a6,t0,t1,cpsii,cpsis,cosp;
                          register __m512 t2,term;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          t1    = _mm512_mul_ps(a,a);
                          cpsii = xcosf(psii);
                          cpsii = _mm512_mul_ps(hlf,cpsii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          cpsis = xcosf(psis);
                          term  = _mm512_fmadd_ps(cpsis,cpsii,cosp);
                          term  = _mm512_mul_ps(term,term);
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,term));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4325_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512  k0  = _mm512_load_ps(&pk0[0]);
                          register __m512  a   = _mm512_load_ps(&pa[0]);
                          register __m512  psis= _mm512_load_ps(&ppsis[0]);
                          register __m512  psii= _mm512_load_ps(&ppsii[0]);
                          register __m512  phi = _mm512_load_ps(&pphi[0]);
                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          const __m512 hlf    = _mm512_set1_ps(0.5f);
                          register __m512 rcs,k04,a6,t0,t1,cpsii,cpsis,cosp;
                          register __m512 t2,term;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          t1    = _mm512_mul_ps(a,a);
                          cpsii = xcosf(psii);
                          cpsii = _mm512_mul_ps(hlf,cpsii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          cpsis = xcosf(psis);
                          term  = _mm512_fmadd_ps(cpsis,cpsii,cosp);
                          term  = _mm512_mul_ps(term,term);
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,term));
                          return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4325_zmm16r4_u(  const float * __restrict  pk0,
                                                const float * __restrict  pa,
                                                const float * __restrict  ppsis,
                                                const float * __restrict  ppsii,
                                                const float * __restrict  pphi) {

                          register __m512  k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512  a   = _mm512_loadu_ps(&pa[0]);
                          register __m512  psis= _mm512_loadu_ps(&ppsis[0]);
                          register __m512  psii= _mm512_loadu_ps(&ppsii[0]);
                          register __m512  phi = _mm512_loadu_ps(&pphi[0]);
                          const __m512 _64pi9 = _mm512_set1_ps(2.263536968418066997601902412409f);
                          const __m512 hlf    = _mm512_set1_ps(0.5f);
                          register __m512 rcs,k04,a6,t0,t1,cpsii,cpsis,cosp;
                          register __m512 t2,term;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_ps(k0,k0);
                          t1    = _mm512_mul_ps(a,a);
                          cpsii = xcosf(psii);
                          cpsii = _mm512_mul_ps(hlf,cpsii);
                          k04   = _mm512_mul_ps(t0,t0);
                          t2    = _mm512_mul_ps(_64pi9,_mm512_mul_ps(k04,k04));
                          a6    = _mm512_mul_ps(t1,_mm512_mul_ps(t1,t1));
                          cpsis = xcosf(psis);
                          term  = _mm512_fmadd_ps(cpsis,cpsii,cosp);
                          term  = _mm512_mul_ps(term,term);
                          rcs   = _mm512_mul_ps(t2,_mm512_mul_ps(a6,term));
                          return (rcs);
                }


                   /*
                          Backscattering RCS for perfectly conducting wire.
                          (2*h>gamma/4)
                          Formula 4.3-29

                     */

                     /*
                          Parameter a1,a2,a3 of equation 4.3-29
                          Formula 4.3-30
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a1_f4330_zmm16r4(const __m512 k0h,
                                           const __m512 psi) {

                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 a1,_2k0h,spsi,arg,_2spsi,sarg;
                          _2k0h = _mm512_add_ps(k0h,k0h);
                          spsi  = xsinf(psi);
                          arg  = _mm512_mul_ps(_2k0h,spsi);
                          _2spsi = _mm512_add_ps(spsi,spsi); 
                          sarg   = xsinf(arg);
                          a1 = _mm512_div_ps(sarg,_2spsi);
                          return (a1); 
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a1_f4330_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 k0h = _mm512_load_ps(&pk0h[0]);
                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 a1,_2k0h,spsi,arg,_2spsi,sarg;
                          _2k0h = _mm512_add_ps(k0h,k0h);
                          spsi  = xsinf(psi);
                          arg  = _mm512_mul_ps(_2k0h,spsi);
                          _2spsi = _mm512_add_ps(spsi,spsi); 
                          sarg   = xsinf(arg);
                          a1 = _mm512_div_ps(sarg,_2spsi);
                          return (a1); 
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a1_f4330_zmm16r4_u(const float * __restrict  pk0h,
                                             const float * __restrict  ppsi) {

                          register __m512 k0h = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 a1,_2k0h,spsi,arg,_2spsi,sarg;
                          _2k0h = _mm512_add_ps(k0h,k0h);
                          spsi  = xsinf(psi);
                          arg  = _mm512_mul_ps(_2k0h,spsi);
                          _2spsi = _mm512_add_ps(spsi,spsi); 
                          sarg   = xsinf(arg);
                          a1 = _mm512_div_ps(sarg,_2spsi);
                          return (a1); 
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a2_f4330_zmm16r4(const __m512 k0h,
                                           const __m512 psi) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_sub_ps(_1,spsi);
                          arg  = _mm512_mul_ps(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_ps(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a2_f4330_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const float * __restrict __ATTR_ALIGN__(64) ppsi) {
   
                          register __m512 k0h = _mm512_load_ps(&pk0h[0]);
                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_sub_ps(_1,spsi);
                          arg  = _mm512_mul_ps(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_ps(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a2_f4330_zmm16r4_u(const float * __restrict  pk0h,
                                             const float * __restrict  ppsi) {
   
                          register __m512 k0h = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_sub_ps(_1,spsi);
                          arg  = _mm512_mul_ps(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_ps(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a3_f4330_zmm16r4(const __m512 k0h,
                                           const __m512 psi) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_add_ps(_1,spsi);
                          arg  = _mm512_mul_ps(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_ps(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a3_f4330_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const float * __restrict __ATTR_ALIGN__(64) ppsi) {
   
                          register __m512 k0h = _mm512_load_ps(&pk0h[0]);
                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_add_ps(_1,spsi);
                          arg  = _mm512_mul_ps(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_ps(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 a3_f4330_zmm16r4_u(const float * __restrict  pk0h,
                                             const float * __restrict  ppsi) {
   
                          register __m512 k0h = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_add_ps(_1,spsi);
                          arg  = _mm512_mul_ps(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_ps(sarg,_1map);
                          return (a2);
                }


                    /*
                          Parameter F1,F2 of equation 4.3-29
                          Formula 4.3-31
                      */

                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F1_f4331_zmm16r4(const __m512 k0a) {

                          const __m512 c0 = _mm512_set1_ps(0.8905f);
                          const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          register __m512 F1,om,om2,arg,larg;
                          arg = _mm512_mul_ps(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_ps(n2,larg);
                          om2 = _mm512_mul_ps(om,om);
                          F1  = _mm512_div_ps(om,_mm512_add_ps(om2,spi));
                          return (F1);
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F1_f4331_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          const __m512 c0 = _mm512_set1_ps(0.8905f);
                          const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          register __m512 F1,om,om2,arg,larg;
                          arg = _mm512_mul_ps(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_ps(n2,larg);
                          om2 = _mm512_mul_ps(om,om);
                          F1  = _mm512_div_ps(om,_mm512_add_ps(om2,spi));
                          return (F1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F1_f4331_zmm16r4_u(const float * __restrict  pk0a) {

                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 c0 = _mm512_set1_ps(0.8905f);
                          const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          register __m512 F1,om,om2,arg,larg;
                          arg = _mm512_mul_ps(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_ps(n2,larg);
                          om2 = _mm512_mul_ps(om,om);
                          F1  = _mm512_div_ps(om,_mm512_add_ps(om2,spi));
                          return (F1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F2_f4331_zmm16r4(const __m512 k0a) {

                          const __m512 c0 = _mm512_set1_ps(0.8905f);
                          const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          register __m512 F1,om,om2,arg,larg;
                          arg = _mm512_mul_ps(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_ps(n2,larg);
                          om2 = _mm512_mul_ps(om,om);
                          F1  = _mm512_div_ps(PI,_mm512_add_ps(om2,spi));
                          return (F1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F2_f4331_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          const __m512 c0 = _mm512_set1_ps(0.8905f);
                          const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          register __m512 F1,om,om2,arg,larg;
                          arg = _mm512_mul_ps(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_ps(n2,larg);
                          om2 = _mm512_mul_ps(om,om);
                          F1  = _mm512_div_ps(PI,_mm512_add_ps(om2,spi));
                          return (F1);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F2_f4331_zmm16r4_u(const float * __restrict  pk0a) {

                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 c0 = _mm512_set1_ps(0.8905f);
                          const __m512 spi= _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 n2 = _mm512_set1_ps(-2.0f);
                          register __m512 F1,om,om2,arg,larg;
                          arg = _mm512_mul_ps(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_ps(n2,larg);
                          om2 = _mm512_mul_ps(om,om);
                          F1  = _mm512_div_ps(PI,_mm512_add_ps(om2,spi));
                          return (F1);
                }


                     /*
                          Parameter (helper) Lambda of equation 4.3-29
                          Formula 4.3-34
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 L_f4334_zmm16r4(const __m512 k0h,
                                          const __m512 k0a) {

                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 c0  = _mm512_set1_ps(4.11f);
                          const __m512 hlf = _mm512_set1_ps(-0.5f);
                          const __m512 n2  = _mm512_set1_ps(-2.0f);
                          const __m512 c1  = _mm512_set1_ps(0.8905f);
                          const __m512 _0  = _mm512_setzero_ps();
                          register __m512 L,om,del,ck0h,sk0h,t0;
                          register __m512 ar1,ar2,lar1,lar2;
                          ar1  = _mm512_mul_ps(k0a,c1);
                          lar1 =xlogf(ar1);
                          ar2  = _mm512_div_ps(k0h,c0);
                          lar2 = xlogf(ar2);
                          om   = _mm512_mul_ps(n2,lar1);
                          del  = _mm512_mul_ps(hlf,lar2);
                          ck0h = xcosf(k0h);
                          t0   = _mm512_sub_ps(_0,_mm512sub_ps(om,del));
                          sk0h = xsinf(k0h);
                          L    = _mm512_fmadd_ps(pi4,sk0h,_mm512_mul_ps(ck0h,t0));
                          return (L);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 L_f4334_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512 k0h = _mm512_load_ps(&pk0h[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 c0  = _mm512_set1_ps(4.11f);
                          const __m512 hlf = _mm512_set1_ps(-0.5f);
                          const __m512 n2  = _mm512_set1_ps(-2.0f);
                          const __m512 c1  = _mm512_set1_ps(0.8905f);
                          const __m512 _0  = _mm512_setzero_ps();
                          register __m512 L,om,del,ck0h,sk0h,t0;
                          register __m512 ar1,ar2,lar1,lar2;
                          ar1  = _mm512_mul_ps(k0a,c1);
                          lar1 =xlogf(ar1);
                          ar2  = _mm512_div_ps(k0h,c0);
                          lar2 = xlogf(ar2);
                          om   = _mm512_mul_ps(n2,lar1);
                          del  = _mm512_mul_ps(hlf,lar2);
                          ck0h = xcosf(k0h);
                          t0   = _mm512_sub_ps(_0,_mm512sub_ps(om,del));
                          sk0h = xsinf(k0h);
                          L    = _mm512_fmadd_ps(pi4,sk0h,_mm512_mul_ps(ck0h,t0));
                          return (L);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 L_f4334_zmm16r4_u(const float * __restrict pk0h,
                                            const float * __restrict pk0a) {

                          register __m512 k0h = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 c0  = _mm512_set1_ps(4.11f);
                          const __m512 hlf = _mm512_set1_ps(-0.5f);
                          const __m512 n2  = _mm512_set1_ps(-2.0f);
                          const __m512 c1  = _mm512_set1_ps(0.8905f);
                          const __m512 _0  = _mm512_setzero_ps();
                          register __m512 L,om,del,ck0h,sk0h,t0;
                          register __m512 ar1,ar2,lar1,lar2;
                          ar1  = _mm512_mul_ps(k0a,c1);
                          lar1 =xlogf(ar1);
                          ar2  = _mm512_div_ps(k0h,c0);
                          lar2 = xlogf(ar2);
                          om   = _mm512_mul_ps(n2,lar1);
                          del  = _mm512_mul_ps(hlf,lar2);
                          ck0h = xcosf(k0h);
                          t0   = _mm512_sub_ps(_0,_mm512sub_ps(om,del));
                          sk0h = xsinf(k0h);
                          L    = _mm512_fmadd_ps(pi4,sk0h,_mm512_mul_ps(ck0h,t0));
                          return (L);
                }


                   /*
                          Parameter (helper) Sigma of equation 4.3-29
                          Formula 4.3-35
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 S_f4335_zmm16r4(const __m512 k0a,
                                          const __m512 k0h) {

                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(7.12f);
                          register __m512 ar,lar,sk0h,ck0h;
                          register __m512 S,t0;
                          ar  = _mm512_mul_ps(c0,k0a);
                          lar = xlogf(ar);
                          sk0h= xsinf(k0h);
                          ck0h= xcosf(k0h);
                          t0  = _mm512_mul_ps(hlf,lar);
                          S   = _mm512_fmsub_ps(t0,sk0h,
                                            _mm512_mul_ps(pi4,ck0h));
                          return (S);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 S_f4335_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0h) {

                          register __m512 k0h = _mm512_load_ps(&pk0h[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(7.12f);
                          register __m512 ar,lar,sk0h,ck0h;
                          register __m512 S,t0;
                          ar  = _mm512_mul_ps(c0,k0a);
                          lar = xlogf(ar);
                          sk0h= xsinf(k0h);
                          ck0h= xcosf(k0h);
                          t0  = _mm512_mul_ps(hlf,lar);
                          S   = _mm512_fmsub_ps(t0,sk0h,
                                            _mm512_mul_ps(pi4,ck0h));
                          return (S);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 S_f4335_zmm16r4_u(const float * __restrict  pk0a,
                                            const float * __restrict  pk0h) {

                          register __m512 k0h = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(7.12f);
                          register __m512 ar,lar,sk0h,ck0h;
                          register __m512 S,t0;
                          ar  = _mm512_mul_ps(c0,k0a);
                          lar = xlogf(ar);
                          sk0h= xsinf(k0h);
                          ck0h= xcosf(k0h);
                          t0  = _mm512_mul_ps(hlf,lar);
                          S   = _mm512_fmsub_ps(t0,sk0h,
                                            _mm512_mul_ps(pi4,ck0h));
                          return (S);
                }


                  /*

                           Parameter G1,G2 of equation 4.3-29
                           Formula 4.3-32
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G2_f4332_zmm16r4(const __m512 k0h,
                                           const __m512 k0a) {

                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 G2,L,S,num,den;
                          L = L_f4334_zmm16r4(k0h,k0a);
                          S = S_f4335_zmm16r4(k0a,k0h);
                          num = _mm512_mul_ps(hlf,S);
                          den = _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          G2  = _mm512_div_ps(num,den);
                          return (G2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G2_f4332_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 G2,L,S,num,den;
                          L = L_f4334_zmm16r4_a(pk0h,pk0a);
                          S = S_f4335_zmm16r4_a(pk0a,pk0h);
                          num = _mm512_mul_ps(hlf,S);
                          den = _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          G2  = _mm512_div_ps(num,den);
                          return (G2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G2_f4332_zmm16r4_u(const float * __restrict  pk0h,
                                             const float * __restrict  pk0a) {

                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 G2,L,S,num,den;
                          L = L_f4334_zmm16r4_a(pk0h,pk0a);
                          S = S_f4335_zmm16r4_a(pk0a,pk0h);
                          num = _mm512_mul_ps(hlf,S);
                          den = _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          G2  = _mm512_div_ps(num,den);
                          return (G2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G1_f4332_zmm16r4(const __m512 k0h,
                                           const __m512 k0a) {

                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 n2  = _mm512_set1_ps(-2.0f);
                         const __m512 c0  = _mm512_set1_ps(0.8905f);
                         register __m512 G1,L,S,om,G1,ln,num,den,om2,t0,rat;
                         L = L_f4334_zmm16r4(k0h,k0a);
                         S = S_f4335_zmm16r4(k0a,k0h);
                         ln= xlogf(_mm512_mul_ps(k0a,c0));
                         om= _mm512_mul_ps(m2,ln); 
                         G2= G2_f4332_zmm16r4(k0h,k0a);
                         om2= _mm512_add_ps(om,om);
                         num= _mm512_mul_ps(hlf,L);
                         t0 = _mm512_mul_ps(PI,G2);
                         ln = _mm512_div_ps(t0,om2);
                         den= _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                         rat= _mm512_div_ps(num,den);
                         G1 = _mm512_sub_ps(rat,ln);
                         return (G1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G1_f4332_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                         register __m512 k0h  = _mm512_load_ps(&pk0h[0]);
                         register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 n2  = _mm512_set1_ps(-2.0f);
                         const __m512 c0  = _mm512_set1_ps(0.8905f);
                         register __m512 G1,L,S,om,G1,ln,num,den,om2,t0,rat;
                         L = L_f4334_zmm16r4(k0h,k0a);
                         S = S_f4335_zmm16r4(k0a,k0h);
                         ln= xlogf(_mm512_mul_ps(k0a,c0));
                         om= _mm512_mul_ps(m2,ln); 
                         G2= G2_f4332_zmm16r4(k0h,k0a);
                         om2= _mm512_add_ps(om,om);
                         num= _mm512_mul_ps(hlf,L);
                         t0 = _mm512_mul_ps(PI,G2);
                         ln = _mm512_div_ps(t0,om2);
                         den= _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                         rat= _mm512_div_ps(num,den);
                         G1 = _mm512_sub_ps(rat,ln);
                         return (G1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G1_f4332_zmm16r4_u(const float * __restrict  pk0h,
                                             const float * __restrict  pk0a) {

                         register __m512 k0h  = _mm512_loadu_ps(&pk0h[0]);
                         register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                         const __m512 hlf = _mm512_set1_ps(0.5f);
                         const __m512 n2  = _mm512_set1_ps(-2.0f);
                         const __m512 c0  = _mm512_set1_ps(0.8905f);
                         register __m512 G1,L,S,om,G1,ln,num,den,om2,t0,rat;
                         L = L_f4334_zmm16r4(k0h,k0a);
                         S = S_f4335_zmm16r4(k0a,k0h);
                         ln= xlogf(_mm512_mul_ps(k0a,c0));
                         om= _mm512_mul_ps(m2,ln); 
                         G2= G2_f4332_zmm16r4(k0h,k0a);
                         om2= _mm512_add_ps(om,om);
                         num= _mm512_mul_ps(hlf,L);
                         t0 = _mm512_mul_ps(PI,G2);
                         ln = _mm512_div_ps(t0,om2);
                         den= _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                         rat= _mm512_div_ps(num,den);
                         G1 = _mm512_sub_ps(rat,ln);
                         return (G1);
                 }


                     /*

                           Parameter H1,H2 of equation 4.3-29
                           Formula 4.3-33
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 H2_f4333_zmm16r4(const __m512 k0h,
                                           const __m512 k0a) {

                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 H2,L,S,num,den,arg;
                          arg  = _mm512_mul_ps(pi2,k0h);
                          L    = L_f4334_zmm16r4(k0h,k0a);
                          S    = S_f4335_zmm16r4(k0a,k0h);
                          num  = _mm512_mul_ps(hlf,S);
                          den  = _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          H2  = _mm512_div_ps(num,den);
                          return (H2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 H2_f4333_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512 k0h  = _mm512_load_ps(&pk0h[0]);
                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 H2,L,S,num,den,arg;
                          arg  = _mm512_mul_ps(pi2,k0h);
                          L    = L_f4334_zmm16r4(k0h,k0a);
                          S    = S_f4335_zmm16r4(k0a,k0h);
                          num  = _mm512_mul_ps(hlf,S);
                          den  = _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          H2  = _mm512_div_ps(num,den);
                          return (H2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 H2_f4333_zmm16r4_u(const float * __restrict  pk0h,
                                             const float * __restrict  pk0a) {

                          register __m512 k0h  = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 H2,L,S,num,den,arg;
                          arg  = _mm512_mul_ps(pi2,k0h);
                          L    = L_f4334_zmm16r4(k0h,k0a);
                          S    = S_f4335_zmm16r4(k0a,k0h);
                          num  = _mm512_mul_ps(hlf,S);
                          den  = _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          H2  = _mm512_div_ps(num,den);
                          return (H2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 H1_f4333_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const float * __restrict __ATTR_ALIGN__(64) pk0h) {

                          register __m512 k0h  = _mm512_load_ps(&pk0h[0]);
                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 n2  = _mm512_set1_ps(-2.0f);
                          const __m512 c0  = _mm512_set1_ps(0.8905f);
                          register __m512 H1,H2,om,ar,lar,L,S,num,den;
                          register __m512 om2,t0,arg;
                          ar = _mm512_mul_ps(k0a,c0);
                          arg= _mm512_mul_ps(k0h,pi2);
                          lar= xlogf(ar);
                          om = _mm512_mul_ps(n2,lar);
                          L  = L_f4334_zmm16r4(k0h,k0a);
                          om2= _mm512_add_ps(om,om);
                          S  = S_f4335_zmm16r4(k0a,k0h);
                          H1 = H1_f4333_zmm16r4(k0h,k0a);
                          num= _mm512_mul_ps(hlf,L);
                          t0 = _mm512_div_ps(_mm512_mul_ps(PI,H1),om2);
                          den= _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          ar = _mm512_div_ps(num,den);
                          H1 = _mm512_sub_ps(ar,t0);
                          return (H1);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 H1_f4333_zmm16r4_u(const float * __restrict  pk0a,
                                             const float * __restrict  pk0h) {

                          register __m512 k0h  = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 n2  = _mm512_set1_ps(-2.0f);
                          const __m512 c0  = _mm512_set1_ps(0.8905f);
                          register __m512 H1,H2,om,ar,lar,L,S,num,den;
                          register __m512 om2,t0,arg;
                          ar = _mm512_mul_ps(k0a,c0);
                          arg= _mm512_mul_ps(k0h,pi2);
                          lar= xlogf(ar);
                          om = _mm512_mul_ps(n2,lar);
                          L  = L_f4334_zmm16r4(k0h,k0a);
                          om2= _mm512_add_ps(om,om);
                          S  = S_f4335_zmm16r4(k0a,k0h);
                          H1 = H1_f4333_zmm16r4(k0h,k0a);
                          num= _mm512_mul_ps(hlf,L);
                          t0 = _mm512_div_ps(_mm512_mul_ps(PI,H1),om2);
                          den= _mm512_fmadd_ps(L,L,_mm512_mul_ps(S,S));
                          ar = _mm512_div_ps(num,den);
                          H1 = _mm512_sub_ps(ar,t0);
                          return (H1);
               }


                 /*
                          Backscattering RCS for perfectly conducting wire.
                          (2*h>gamma/4)
                          Formula 4.3-29

                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4329_zmm16r4(const __m512 k0,
                                            const __m512 gami,
                                            const __m512 gams,
                                            const __m512 k0h,
                                            const __m512 k0a,
                                            const __m512 psi) {

                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          const __m512 _2    = _mm512_set1_ps(2.0f);
                          register __m512 rcs,a1,a2,a3,F1,F2,G1,G2,H1,H2,first;
                          register __m512 cgami,cgams,c2gami,c2gams,sinps;
                          register __m512 arg,sarg,carg,t0,t1,t2,t3,t4,x0,x1,t5,b0;
                          register __m512 a1s,F1F2,G1G2,a2pa3,a2ma3,H1H2,a2sma3s;
                          register __m512 GHGH,_2a1,FGFG,FHFH,tmp1,tmp2,tmp3;
                          b0     = _mm512_div_ps(_16pi,_mm512_mul_ps(k0,k0));
                          a1     = a1_f4330_zmm16r4(k0h,psi);
                          _2a1   = _mm512_add_ps(a1,a1);
                          cgami  = xcosf(gami);
                          F1     = F1_f4331_zmm16r4(k0a);
                          c2gami = _mm512_mul_ps(cgami,cgami);
                          F2     = F2_f4331_zmm16r4(k0a);
                          cgams  = xcosf(gams);
                          G1     = G1_f4332_zmm16r4(k0h,k0a);
                          c2gams = _mm512_mul_ps(cgams,cgams);
                          a2     = a2_f4330_zmm16r4(k0h,psi);
                          first  = _mm512_mul_ps(b0,_mm512_mul_ps(c2gami,c2gams));
                          G2     = G1_f4332_zmm16r4(k0h,k0a);
                          sinps  = xsinf(psi);
                          a3     = a3_f4330_zmm16r4(k0h,psi);
                          H1     = H1_f4333_zmm16r4(k0h,k0a);
                          arg    = _mm512_mul_ps(k0h,sinps);
                          H2     = H2_f4333_zmm16r4(k0h,k0a);
                          sarg   = xsinf(arg);
                          a1s    = _mm512_mul_ps(a1,a1);
                          carg   = xcosf(arg);
                          x0     = _mm512_add_ps(a2,a3);
                          a2pa3  = _mm512_mul_ps(x0,x0);
                          F1F2   = _mm512_fmadd_ps(F1,F1,_mm512_mul_ps(F2,F2));
                          x1     = _mm512_sub_ps(a2,a3);
                          t0     = _mm512_mul_ps(a1s,F1F2);
                          a2ma3  = _mm512_mul_ps(x1,x1);
                          G1G2   = _mm512_fmadd_ps(G1,G1,_mm512_mul_ps(G2,G2));
                          t1     = _mm512_mul_ps(a2pa3,_mm512_mul_ps(G1G2,carg));
                          x0     = _mm512_mul_ps(sarg,sarg);
                          H1H2   = _mm512_fmadd_ps(H1,H1,_mm512_mul_ps(H2,H2));
                          t2     = _mm512_mul_ps(a2ma3,_mm512_mul_ps(H1H2,x0));
                          a2sma3s= _mm512_mul_ps(_2,_mm512_fmsub_ps(a2,a2,
                                                                _mm512_mul_ps(a3,a3)));
                          GHGH   = _mm512_fmadd_ps(G1,H1,_mm512_mul_ps(G2,H2));
                          x1     = _mm512_mul_ps(carg,sarg);
                          t3     = _mm512_mul_ps(a2sma3s,_mm512_mul_ps(GHGH,x1));
                          x0     = _mm512_mul_ps(_2a1,a2pa3);
                          FGFG   = _mm512_fmadd_ps(F1,G1,_mm512_mul_ps(F2,G2));
                          t4     = _mm512_mul_ps(x0,_mm512_mul_ps(FGFG,carg);
                          x1     = _mm512_mul_ps(_2a1,a2ma3);
                          FHFH   = _mm512_fmadd_ps(F1,H1,_mm512_mul_ps(F2,H2));
                          t5     = _mm512_mul_ps(x1,_mm512_mul_ps(FHFH,sarg));
                          tmp1   = _mm512_add_ps(t0,_mm512_add_ps(t1,t2));
                          tmp2   = _mm512_sub_ps(_mm512_add_ps(t3,t4),t5);
                          tmp3   = _mm512_sub_ps(tmp1,tmp2);
                          rcs    = _mm512_mul_ps(first,tmp3);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4329_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pgami,
                                              const float * __restrict __ATTR_ALIGN__(64) pgams,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0h,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsi) {


                          register __m512 k0     = _mm512_load_ps(&pk0[0]);
                          register __m512 gami   = _mm512_load_ps(&pgami[0]);
                          register __m512 gams   = _mm512_load_ps(&pgams[0]);
                          register __m512 k0h    = _mm512_load_ps(&pk0h[0]);
                          register __m512 k0a    = _mm512_load_ps(&pk0a[0]);
                          register __m512 psi    = _mm512_load_ps(&ppsi[0]);  
                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          const __m512 _2    = _mm512_set1_ps(2.0f);
                          register __m512 rcs,a1,a2,a3,F1,F2,G1,G2,H1,H2,first;
                          register __m512 cgami,cgams,c2gami,c2gams,sinps;
                          register __m512 arg,sarg,carg,t0,t1,t2,t3,t4,x0,x1,t5,b0;
                          register __m512 a1s,F1F2,G1G2,a2pa3,a2ma3,H1H2,a2sma3s;
                          register __m512 GHGH,_2a1,FGFG,FHFH,tmp1,tmp2,tmp3;
                          b0     = _mm512_div_ps(_16pi,_mm512_mul_ps(k0,k0));
                          a1     = a1_f4330_zmm16r4(k0h,psi);
                          _2a1   = _mm512_add_ps(a1,a1);
                          cgami  = xcosf(gami);
                          F1     = F1_f4331_zmm16r4(k0a);
                          c2gami = _mm512_mul_ps(cgami,cgami);
                          F2     = F2_f4331_zmm16r4(k0a);
                          cgams  = xcosf(gams);
                          G1     = G1_f4332_zmm16r4(k0h,k0a);
                          c2gams = _mm512_mul_ps(cgams,cgams);
                          a2     = a2_f4330_zmm16r4(k0h,psi);
                          first  = _mm512_mul_ps(b0,_mm512_mul_ps(c2gami,c2gams));
                          G2     = G1_f4332_zmm16r4(k0h,k0a);
                          sinps  = xsinf(psi);
                          a3     = a3_f4330_zmm16r4(k0h,psi);
                          H1     = H1_f4333_zmm16r4(k0h,k0a);
                          arg    = _mm512_mul_ps(k0h,sinps);
                          H2     = H2_f4333_zmm16r4(k0h,k0a);
                          sarg   = xsinf(arg);
                          a1s    = _mm512_mul_ps(a1,a1);
                          carg   = xcosf(arg);
                          x0     = _mm512_add_ps(a2,a3);
                          a2pa3  = _mm512_mul_ps(x0,x0);
                          F1F2   = _mm512_fmadd_ps(F1,F1,_mm512_mul_ps(F2,F2));
                          x1     = _mm512_sub_ps(a2,a3);
                          t0     = _mm512_mul_ps(a1s,F1F2);
                          a2ma3  = _mm512_mul_ps(x1,x1);
                          G1G2   = _mm512_fmadd_ps(G1,G1,_mm512_mul_ps(G2,G2));
                          t1     = _mm512_mul_ps(a2pa3,_mm512_mul_ps(G1G2,carg));
                          x0     = _mm512_mul_ps(sarg,sarg);
                          H1H2   = _mm512_fmadd_ps(H1,H1,_mm512_mul_ps(H2,H2));
                          t2     = _mm512_mul_ps(a2ma3,_mm512_mul_ps(H1H2,x0));
                          a2sma3s= _mm512_mul_ps(_2,_mm512_fmsub_ps(a2,a2,
                                                                _mm512_mul_ps(a3,a3)));
                          GHGH   = _mm512_fmadd_ps(G1,H1,_mm512_mul_ps(G2,H2));
                          x1     = _mm512_mul_ps(carg,sarg);
                          t3     = _mm512_mul_ps(a2sma3s,_mm512_mul_ps(GHGH,x1));
                          x0     = _mm512_mul_ps(_2a1,a2pa3);
                          FGFG   = _mm512_fmadd_ps(F1,G1,_mm512_mul_ps(F2,G2));
                          t4     = _mm512_mul_ps(x0,_mm512_mul_ps(FGFG,carg);
                          x1     = _mm512_mul_ps(_2a1,a2ma3);
                          FHFH   = _mm512_fmadd_ps(F1,H1,_mm512_mul_ps(F2,H2));
                          t5     = _mm512_mul_ps(x1,_mm512_mul_ps(FHFH,sarg));
                          tmp1   = _mm512_add_ps(t0,_mm512_add_ps(t1,t2));
                          tmp2   = _mm512_sub_ps(_mm512_add_ps(t3,t4),t5);
                          tmp3   = _mm512_sub_ps(tmp1,tmp2);
                          rcs    = _mm512_mul_ps(first,tmp3);
                          return (rcs);
               }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4329_zmm16r4_u(const float * __restrict pk0,
                                              const float * __restrict  pgami,
                                              const float * __restrict  pgams,
                                              const float * __restrict  pk0h,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  ppsi) {


                          register __m512 k0     = _mm512_loadu_ps(&pk0[0]);
                          register __m512 gami   = _mm512_loadu_ps(&pgami[0]);
                          register __m512 gams   = _mm512_loadu_ps(&pgams[0]);
                          register __m512 k0h    = _mm512_loadu_ps(&pk0h[0]);
                          register __m512 k0a    = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 psi    = _mm512_loadu_ps(&ppsi[0]);  
                          const __m512 _16pi = _mm512_set1_ps(50.265482457436691815402294132472f);
                          const __m512 _2    = _mm512_set1_ps(2.0f);
                          register __m512 rcs,a1,a2,a3,F1,F2,G1,G2,H1,H2,first;
                          register __m512 cgami,cgams,c2gami,c2gams,sinps;
                          register __m512 arg,sarg,carg,t0,t1,t2,t3,t4,x0,x1,t5,b0;
                          register __m512 a1s,F1F2,G1G2,a2pa3,a2ma3,H1H2,a2sma3s;
                          register __m512 GHGH,_2a1,FGFG,FHFH,tmp1,tmp2,tmp3;
                          b0     = _mm512_div_ps(_16pi,_mm512_mul_ps(k0,k0));
                          a1     = a1_f4330_zmm16r4(k0h,psi);
                          _2a1   = _mm512_add_ps(a1,a1);
                          cgami  = xcosf(gami);
                          F1     = F1_f4331_zmm16r4(k0a);
                          c2gami = _mm512_mul_ps(cgami,cgami);
                          F2     = F2_f4331_zmm16r4(k0a);
                          cgams  = xcosf(gams);
                          G1     = G1_f4332_zmm16r4(k0h,k0a);
                          c2gams = _mm512_mul_ps(cgams,cgams);
                          a2     = a2_f4330_zmm16r4(k0h,psi);
                          first  = _mm512_mul_ps(b0,_mm512_mul_ps(c2gami,c2gams));
                          G2     = G1_f4332_zmm16r4(k0h,k0a);
                          sinps  = xsinf(psi);
                          a3     = a3_f4330_zmm16r4(k0h,psi);
                          H1     = H1_f4333_zmm16r4(k0h,k0a);
                          arg    = _mm512_mul_ps(k0h,sinps);
                          H2     = H2_f4333_zmm16r4(k0h,k0a);
                          sarg   = xsinf(arg);
                          a1s    = _mm512_mul_ps(a1,a1);
                          carg   = xcosf(arg);
                          x0     = _mm512_add_ps(a2,a3);
                          a2pa3  = _mm512_mul_ps(x0,x0);
                          F1F2   = _mm512_fmadd_ps(F1,F1,_mm512_mul_ps(F2,F2));
                          x1     = _mm512_sub_ps(a2,a3);
                          t0     = _mm512_mul_ps(a1s,F1F2);
                          a2ma3  = _mm512_mul_ps(x1,x1);
                          G1G2   = _mm512_fmadd_ps(G1,G1,_mm512_mul_ps(G2,G2));
                          t1     = _mm512_mul_ps(a2pa3,_mm512_mul_ps(G1G2,carg));
                          x0     = _mm512_mul_ps(sarg,sarg);
                          H1H2   = _mm512_fmadd_ps(H1,H1,_mm512_mul_ps(H2,H2));
                          t2     = _mm512_mul_ps(a2ma3,_mm512_mul_ps(H1H2,x0));
                          a2sma3s= _mm512_mul_ps(_2,_mm512_fmsub_ps(a2,a2,
                                                                _mm512_mul_ps(a3,a3)));
                          GHGH   = _mm512_fmadd_ps(G1,H1,_mm512_mul_ps(G2,H2));
                          x1     = _mm512_mul_ps(carg,sarg);
                          t3     = _mm512_mul_ps(a2sma3s,_mm512_mul_ps(GHGH,x1));
                          x0     = _mm512_mul_ps(_2a1,a2pa3);
                          FGFG   = _mm512_fmadd_ps(F1,G1,_mm512_mul_ps(F2,G2));
                          t4     = _mm512_mul_ps(x0,_mm512_mul_ps(FGFG,carg);
                          x1     = _mm512_mul_ps(_2a1,a2ma3);
                          FHFH   = _mm512_fmadd_ps(F1,H1,_mm512_mul_ps(F2,H2));
                          t5     = _mm512_mul_ps(x1,_mm512_mul_ps(FHFH,sarg));
                          tmp1   = _mm512_add_ps(t0,_mm512_add_ps(t1,t2));
                          tmp2   = _mm512_sub_ps(_mm512_add_ps(t3,t4),t5);
                          tmp3   = _mm512_sub_ps(tmp1,tmp2);
                          rcs    = _mm512_mul_ps(first,tmp3);
                          return (rcs);
               }


                  /*

                         Simplified back and bistatic scattering RCS for
                         half and full-wave dipole (2*h == gam0/2, and gam0)
                         gam0 -- wavelength.
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4337_zmm16r4(const __m512 gammi,
                                            const __m512 gamms,
                                            const __m512 psii,
                                            const __m512 psis,
                                            const __m512 g0 )  {//wavelength coeff

                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          register __m512 rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512 spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_ps(pi2,spsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_ps(pi2,spsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_ps(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_ps(cgami,cgami);
                          t0    = _mm512_mul_ps(g0,_mm512_mul_ps(c2gami,c2gams));
                          c1    = xcosf(carg1);
                          rat1  = _mm512_div_ps(c1,cpsii);
                          tmp0  = _mm512_mul_ps(rat1,rat1);
                          c2    = xcosf(carg2);
                          rat2  = _mm512_div_ps(c2,cpsis);
                          tmp1  = _mm512_mul_ps(rat2,rat2);
                          t1    = _mm512_mul_ps(tmp0,tmp1);
                          rcs   = _mm512_mul_ps(t0,t1);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4337_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgammi,
                                              const float * __restrict __ATTR_ALIGN__(64) pgamms,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64) pg0 )  {//wavelength coeff

                          register __m512  gammi = _mm512_load_ps(&pgammi[0]);
                          register __m512  gamms = _mm512_load_ps(&pgamms[0]);
                          register __m512  psii  = _mm512_load_ps(&ppsii[0]);
                          register __m512  psis  = _mm512_load_ps(&ppsis[0]);
                          register __m512  g0    = _mm512_load_ps(&pg0[0]);
                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          register __m512 rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512 spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_ps(pi2,spsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_ps(pi2,spsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_ps(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_ps(cgami,cgami);
                          t0    = _mm512_mul_ps(g0,_mm512_mul_ps(c2gami,c2gams));
                          c1    = xcosf(carg1);
                          rat1  = _mm512_div_ps(c1,cpsii);
                          tmp0  = _mm512_mul_ps(rat1,rat1);
                          c2    = xcosf(carg2);
                          rat2  = _mm512_div_ps(c2,cpsis);
                          tmp1  = _mm512_mul_ps(rat2,rat2);
                          t1    = _mm512_mul_ps(tmp0,tmp1);
                          rcs   = _mm512_mul_ps(t0,t1);
                          return (rcs);
                 }
      


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4337_zmm16r4_u(const float * __restrict  pgammi,
                                              const float * __restrict  pgamms,
                                              const float * __restrict  ppsii,
                                              const float * __restrict  ppsis,
                                              const float * __restrict  pg0 )  {//wavelength coeff

                          register __m512  gammi = _mm512_loadu_ps(&pgammi[0]);
                          register __m512  gamms = _mm512_loadu_ps(&pgamms[0]);
                          register __m512  psii  = _mm512_loadu_ps(&ppsii[0]);
                          register __m512  psis  = _mm512_loadu_ps(&ppsis[0]);
                          register __m512  g0    = _mm512_loadu_ps(&pg0[0]);
                          const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                          register __m512 rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512 spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_ps(pi2,spsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_ps(pi2,spsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_ps(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_ps(cgami,cgami);
                          t0    = _mm512_mul_ps(g0,_mm512_mul_ps(c2gami,c2gams));
                          c1    = xcosf(carg1);
                          rat1  = _mm512_div_ps(c1,cpsii);
                          tmp0  = _mm512_mul_ps(rat1,rat1);
                          c2    = xcosf(carg2);
                          rat2  = _mm512_div_ps(c2,cpsis);
                          tmp1  = _mm512_mul_ps(rat2,rat2);
                          t1    = _mm512_mul_ps(tmp0,tmp1);
                          rcs   = _mm512_mul_ps(t0,t1);
                          return (rcs);
                 }


                   /*

                         Simplified back and bistatic scattering RCS for
                         Full-wave dipole (2*h == gam0)
                         gam0 -- wavelength.
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4340_zmm16r4(const __m512 gammi,
                                            const __m512 gamms,
                                            const __m512 psii,
                                            const __m512 psis,
                                            const __m512 g0 )  {//wavelength coeff

                         
                          register __m512 rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512 spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_ps(PI,spsii);
                          cpsii = _mm512_mul_ps(cpsii,cpsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_ps(PI,spsis);
                          cpsis = _mm512_mul_ps(cpsis,cpsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_ps(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_ps(cgami,cgami);
                          t0    = _mm512_mul_ps(g0,_mm512_mul_ps(c2gami,c2gams));
                          c1    = xsinf(carg1);
                          rat1  = _mm512_div_ps(c1,cpsii);
                          tmp0  = _mm512_mul_ps(rat1,rat1);
                          c2    = xsinf(carg2);
                          rat2  = _mm512_div_ps(c2,cpsis);
                          tmp1  = _mm512_mul_ps(rat2,rat2);
                          t1    = _mm512_mul_ps(tmp0,tmp1);
                          rcs   = _mm512_mul_ps(t0,t1);
                          return (rcs);
                 }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4340_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgammi,
                                              const float * __restrict __ATTR_ALIGN__(64) pgamms,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64) pg0 )  {//wavelength coeff

                         
                          register __m512  gammi = _mm512_load_ps(&pgammi[0]);
                          register __m512  gamms = _mm512_load_ps(&pgamms[0]);
                          register __m512  psii  = _mm512_load_ps(&ppsii[0]);
                          register __m512  psis  = _mm512_load_ps(&ppsis[0]);
                          register __m512  g0    = _mm512_load_ps(&pg0[0]);
                          register __m512 rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512 spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_ps(PI,spsii);
                          cpsii = _mm512_mul_ps(cpsii,cpsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_ps(PI,spsis);
                          cpsis = _mm512_mul_ps(cpsis,cpsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_ps(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_ps(cgami,cgami);
                          t0    = _mm512_mul_ps(g0,_mm512_mul_ps(c2gami,c2gams));
                          c1    = xsinf(carg1);
                          rat1  = _mm512_div_ps(c1,cpsii);
                          tmp0  = _mm512_mul_ps(rat1,rat1);
                          c2    = xsinf(carg2);
                          rat2  = _mm512_div_ps(c2,cpsis);
                          tmp1  = _mm512_mul_ps(rat2,rat2);
                          t1    = _mm512_mul_ps(tmp0,tmp1);
                          rcs   = _mm512_mul_ps(t0,t1);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4340_zmm16r4_u(const float * __restrict  pgammi,
                                              const float * __restrict  pgamms,
                                              const float * __restrict  ppsii,
                                              const float * __restrict  ppsis,
                                              const float * __restrict  pg0 )  {//wavelength coeff

                         
                          register __m512  gammi = _mm512_loadu_ps(&pgammi[0]);
                          register __m512  gamms = _mm512_loadu_ps(&pgamms[0]);
                          register __m512  psii  = _mm512_loadu_ps(&ppsii[0]);
                          register __m512  psis  = _mm512_loadu_ps(&ppsis[0]);
                          register __m512  g0    = _mm512_loadu_ps(&pg0[0]);
                          register __m512 rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512 spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_ps(PI,spsii);
                          cpsii = _mm512_mul_ps(cpsii,cpsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_ps(PI,spsis);
                          cpsis = _mm512_mul_ps(cpsis,cpsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_ps(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_ps(cgami,cgami);
                          t0    = _mm512_mul_ps(g0,_mm512_mul_ps(c2gami,c2gams));
                          c1    = xsinf(carg1);
                          rat1  = _mm512_div_ps(c1,cpsii);
                          tmp0  = _mm512_mul_ps(rat1,rat1);
                          c2    = xsinf(carg2);
                          rat2  = _mm512_div_ps(c2,cpsis);
                          tmp1  = _mm512_mul_ps(rat2,rat2);
                          t1    = _mm512_mul_ps(tmp0,tmp1);
                          rcs   = _mm512_mul_ps(t0,t1);
                          return (rcs);
                 }


                     /*
                           Cylinder length much greater then wavelength (h>>gamma).
                           Biscattering RCS, formula 4.3-43
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4343_zmm16r4(const __m512 rcs_inf, // rcs of inifnitely long cylinder (section 4.2)
                                            const __m512 k0,
                                            const __m512 h,
                                            const __m512 psis,
                                            const __m512 psii) {

                           const __m512 _4 = _mm512_set1_ps(4.0f);
                           register __m512 k0h,x0,term1,cpsis,c2psis,rcs;
                           register __m512 term2,spsii,spsis,arg,sarg,rat;
                           k0h   = _mm512_mul_ps(_4,_mm512_mul_ps(k0,h));
                           x0    = _mm512_mul_ps(k0h,k0h);
                           cpsis = xcosf(psis);
                           term1 = _mm512_div_ps(x0,PI);
                           c2psis= _mm512_mul_ps(cpsis,cpsis);
                           term1 = _mm512_mul_ps(term1,_mm512_mul_ps(c2psis,rcs_inf));
                           spsis = xsinf(psis);
                           spsii = xsinf(psii);
                           x0    = _mm512_add_ps(spsis,spsii);
                           arg   = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                           sarg  = xsinf(arg);
                           rat   = _mm512_div_ps(sarg,arg);
                           term2 = _mm512_mul_ps(rat,rat);
                           rcs   = _mm512_mul_ps(term1,term2);
                           return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4343_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) prcs_inf, // rcs of inifnitely long cylinder (section 4.2)
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) ph,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii) {

                           register __m512 prcs_inf  = _mm512_load_ps(&prcs_inf[0]);
                           register __m512 k0        = _mm512_load_ps(&pk0[0]);
                           register __m512 ph        = _mm512_load_ps(&ph[0]);
                           register __m512 psis      = _mm512_load_ps(&ppsis[0]);
                           register __m512 ppsis     = _mm512_load_ps(&ppsii[0]);
                           const __m512 _4 = _mm512_set1_ps(4.0f);
                           register __m512 k0h,x0,term1,cpsis,c2psis,rcs;
                           register __m512 term2,spsii,spsis,arg,sarg,rat;
                           k0h   = _mm512_mul_ps(_4,_mm512_mul_ps(k0,h));
                           x0    = _mm512_mul_ps(k0h,k0h);
                           cpsis = xcosf(psis);
                           term1 = _mm512_div_ps(x0,PI);
                           c2psis= _mm512_mul_ps(cpsis,cpsis);
                           term1 = _mm512_mul_ps(term1,_mm512_mul_ps(c2psis,rcs_inf));
                           spsis = xsinf(psis);
                           spsii = xsinf(psii);
                           x0    = _mm512_add_ps(spsis,spsii);
                           arg   = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                           sarg  = xsinf(arg);
                           rat   = _mm512_div_ps(sarg,arg);
                           term2 = _mm512_mul_ps(rat,rat);
                           rcs   = _mm512_mul_ps(term1,term2);
                           return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4343_zmm16r4_u(const float * __restrict  prcs_inf, // rcs of inifnitely long cylinder (section 4.2)
                                              const float * __restrict  pk0,
                                              const float * __restrict  ph,
                                              const float * __restrict  ppsis,
                                              const float * __restrict  ppsii) {

                           register __m512 prcs_inf  = _mm512_loadu_ps(&prcs_inf[0]);
                           register __m512 k0        = _mm512_loadu_ps(&pk0[0]);
                           register __m512 ph        = _mm512_loadu_ps(&ph[0]);
                           register __m512 psis      = _mm512_loadu_ps(&ppsis[0]);
                           register __m512 ppsis     = _mm512_loadu_ps(&ppsii[0]);
                           const __m512 _4 = _mm512_set1_ps(4.0f);
                           register __m512 k0h,x0,term1,cpsis,c2psis,rcs;
                           register __m512 term2,spsii,spsis,arg,sarg,rat;
                           k0h   = _mm512_mul_ps(_4,_mm512_mul_ps(k0,h));
                           x0    = _mm512_mul_ps(k0h,k0h);
                           cpsis = xcosf(psis);
                           term1 = _mm512_div_ps(x0,PI);
                           c2psis= _mm512_mul_ps(cpsis,cpsis);
                           term1 = _mm512_mul_ps(term1,_mm512_mul_ps(c2psis,rcs_inf));
                           spsis = xsinf(psis);
                           spsii = xsinf(psii);
                           x0    = _mm512_add_ps(spsis,spsii);
                           arg   = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                           sarg  = xsinf(arg);
                           rat   = _mm512_div_ps(sarg,arg);
                           term2 = _mm512_mul_ps(rat,rat);
                           rcs   = _mm512_mul_ps(term1,term2);
                           return (rcs);
                }


                  /*
                         General bistatic scattering RCS from long thin wire.
                         Formula 4.3-44
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4344_zmm16r4(const __m512 h,
                                            const __m512 k0,
                                            const __m512 k0a,
                                            const __m512 psii,
                                            const __m512 psis,
                                            const __m512 gams,
                                            const __m512 gami) {

                          const __m512 c0 = _mm512_set1_ps(12.566370614359172953850573533118f);
                          const __m512 c1 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const __m512 c2 = _mm512_set1_ps(0.8905f);
                          register __m512 term1,term2,term3,cgami,cgams,c2gami,c2gams;
                          register __m512 rcs,inv,arg,sarg,rat1,rat2,x0,x1,arg2,larg;
                          register __m512 cpsii,cpsis,fac,c2psii,c2psis,spsii,spsis;
                          fac    = _mm512_mul_ps(c0,_mm512_mul_ps(h,h));
                          arg2   = _mm512_mul_ps(k0a,c2);
                          cpsii  = xcosf(psii);
                          cpsis  = xcosf(psis);
                          c2psii = _mm512_mul_ps(cpsii,cpsii);
                          c2psis = _mm512_mul_ps(cpsis,cpsis);
                          arg2   = _mm512_mul_ps(cpsii,arg2);
                          rat1   = _mm512_div_ps(c2psis,c2psii);
                          cgami  = xcosf(gami);
                          c2gami = _mm512_mul_ps(cgami,cgami);
                          cgams  = xcosf(gams);
                          c2gams = _mm512_mul_ps(cgams,cgams);
                          x0     = _mm512_mul_ps(c2gams,c2gami);
                          term1  = _mm512_mul_ps(fac,_mm512_mul_ps(rat1,x0));
                          larg   = xlogf(arg2);
                          spsii  = xsinf(psii);
                          x1     = _mm512_fmadd_ps(larg,larg,c1);
                          inv    = _mm512_rcp14_ps(x1);
                          spsis  = xsinf(psis);
                          x0     = _mm512_add_ps(spsii,spsis);
                          arg    = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                          sarg   = xsinf(arg);
                          rat2   = _mm512_div_ps(sarg,arg);
                          term2  = _mm512_mul_ps(rat2,rat2);
                          rcs    = _mm512_mul_ps(term1,_mm512_mul_ps(inv,term2));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4344_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  ph,
                                              const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64)  ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64)  ppsis,
                                              const float * __restrict __ATTR_ALIGN__(64)  pgams,
                                              const float * __restrict __ATTR_ALIGN__(64)  pgami) {

                          register __m512 h    = _mm512_load_ps(&ph[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]); 
                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          register __m512 psii = _mm512_load_ps(&ppsii[0]);
                          register __m512 psis = _mm512_load_ps(&ppsis[0]);
                          register __m512 gams = _mm512_load_ps(&pgams[0]);
                          register __m512 gami = _mm512_load_ps(&pgami[0]);
                          const __m512 c0 = _mm512_set1_ps(12.566370614359172953850573533118f);
                          const __m512 c1 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const __m512 c2 = _mm512_set1_ps(0.8905f);
                          register __m512 term1,term2,term3,cgami,cgams,c2gami,c2gams;
                          register __m512 rcs,inv,arg,sarg,rat1,rat2,x0,x1,arg2,larg;
                          register __m512 cpsii,cpsis,fac,c2psii,c2psis,spsii,spsis;
                          fac    = _mm512_mul_ps(c0,_mm512_mul_ps(h,h));
                          arg2   = _mm512_mul_ps(k0a,c2);
                          cpsii  = xcosf(psii);
                          cpsis  = xcosf(psis);
                          c2psii = _mm512_mul_ps(cpsii,cpsii);
                          c2psis = _mm512_mul_ps(cpsis,cpsis);
                          arg2   = _mm512_mul_ps(cpsii,arg2);
                          rat1   = _mm512_div_ps(c2psis,c2psii);
                          cgami  = xcosf(gami);
                          c2gami = _mm512_mul_ps(cgami,cgami);
                          cgams  = xcosf(gams);
                          c2gams = _mm512_mul_ps(cgams,cgams);
                          x0     = _mm512_mul_ps(c2gams,c2gami);
                          term1  = _mm512_mul_ps(fac,_mm512_mul_ps(rat1,x0));
                          larg   = xlogf(arg2);
                          spsii  = xsinf(psii);
                          x1     = _mm512_fmadd_ps(larg,larg,c1);
                          inv    = _mm512_rcp14_ps(x1);
                          spsis  = xsinf(psis);
                          x0     = _mm512_add_ps(spsii,spsis);
                          arg    = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                          sarg   = xsinf(arg);
                          rat2   = _mm512_div_ps(sarg,arg);
                          term2  = _mm512_mul_ps(rat2,rat2);
                          rcs    = _mm512_mul_ps(term1,_mm512_mul_ps(inv,term2));
                          return (rcs);
                }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4344_zmm16r4_u(const float * __restrict   ph,
                                              const float * __restrict   pk0,
                                              const float * __restrict   pk0a,
                                              const float * __restrict   ppsii,
                                              const float * __restrict   ppsis,
                                              const float * __restrict   pgams,
                                              const float * __restrict   pgami) {

                          register __m512 h    = _mm512_loadu_ps(&ph[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]); 
                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                          register __m512 psis = _mm512_loadu_ps(&ppsis[0]);
                          register __m512 gams = _mm512_loadu_ps(&pgams[0]);
                          register __m512 gami = _mm512_loadu_ps(&pgami[0]);
                          const __m512 c0 = _mm512_set1_ps(12.566370614359172953850573533118f);
                          const __m512 c1 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const __m512 c2 = _mm512_set1_ps(0.8905f);
                          register __m512 term1,term2,term3,cgami,cgams,c2gami,c2gams;
                          register __m512 rcs,inv,arg,sarg,rat1,rat2,x0,x1,arg2,larg;
                          register __m512 cpsii,cpsis,fac,c2psii,c2psis,spsii,spsis;
                          fac    = _mm512_mul_ps(c0,_mm512_mul_ps(h,h));
                          arg2   = _mm512_mul_ps(k0a,c2);
                          cpsii  = xcosf(psii);
                          cpsis  = xcosf(psis);
                          c2psii = _mm512_mul_ps(cpsii,cpsii);
                          c2psis = _mm512_mul_ps(cpsis,cpsis);
                          arg2   = _mm512_mul_ps(cpsii,arg2);
                          rat1   = _mm512_div_ps(c2psis,c2psii);
                          cgami  = xcosf(gami);
                          c2gami = _mm512_mul_ps(cgami,cgami);
                          cgams  = xcosf(gams);
                          c2gams = _mm512_mul_ps(cgams,cgams);
                          x0     = _mm512_mul_ps(c2gams,c2gami);
                          term1  = _mm512_mul_ps(fac,_mm512_mul_ps(rat1,x0));
                          larg   = xlogf(arg2);
                          spsii  = xsinf(psii);
                          x1     = _mm512_fmadd_ps(larg,larg,c1);
                          inv    = _mm512_rcp14_ps(x1);
                          spsis  = xsinf(psis);
                          x0     = _mm512_add_ps(spsii,spsis);
                          arg    = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                          sarg   = xsinf(arg);
                          rat2   = _mm512_div_ps(sarg,arg);
                          term2  = _mm512_mul_ps(rat2,rat2);
                          rcs    = _mm512_mul_ps(term1,_mm512_mul_ps(inv,term2));
                          return (rcs);
                }


                    /*

                          General backscatter (only) scattering RCS from long thin wire.
                          Formula 4.3-45
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4345_zmm16r4(const __m512 psi,
                                            const __m512 k0a,
                                            const __m512 gami,
                                            const __m512 gams,
                                            const __m512 k0,
                                            const __m512 h) {

                         const __m512 pi24 = _mm512_set1_ps(2.467401100272339654708622749969f);
                         const __m512 _2pi = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 c0   = _mm512_set1_ps(0.8905f);
                         register __m512 rat1,arg,sarg,arg2,larg2,k0h,t0,rat;
                         register __m512 rcs,cpsi,cgami,cgams,c2gami,c2gams,spsi;
                         register __m512 x0,x1;
                         k0h   = _mm512_mul_ps(k0,h);
                         t0    = _mm512_mul_ps(_2pi,_mm512_mul_ps(h,h));
                         x0    = _mm512_add_ps(k0h,k0h);
                         spsi  = xsinf(psi);
                         arg   = _mm512_mul_ps(x0,spsi);
                         cpsi  = xcosf(psi);
                         arg2  = _mm512_mul_ps(cpsi,_mm512_mul_ps(k0a,c0));
                         larg  = _mm512_fmadd_ps(arg2,arg2,pi24);
                         sarg  = xsinf(arg);
                         cgams = xcosf(gams);
                         rat   = _mm512_div_ps(sarg,arg);
                         cgami = xcosf(gami);
                         x1    = _mm512_mul_ps(rat,rat);
                         c2gams= _mm512_mul_ps(cgams,cgams);
                         c2gami= _mm512_mul_ps(cgami,cgami);
                         x0    = _mm512_mul_ps(t0,_mm512_mul_ps(c2gams,c2gami));
                         rat1  = _mm512_div_ps(x0,larg);
                         rcs   = _mm512_mul_ps(rat1,x1);
                         return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4345_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  ph,
                                              const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64)  pgams,
                                              const float * __restrict __ATTR_ALIGN__(64)  pgami,
                                              const float * __restrict __ATTR_ALIGN__(64)  ppsi) {

                         register __m512 h    = _mm512_load_ps(&ph[0]);
                         register __m512 k0   = _mm512_load_ps(&pk0[0]); 
                         register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                         register __m512 psi  = _mm512_load_ps(&ppsi[0]);
                         register __m512 gams = _mm512_load_ps(&pgams[0]);
                         register __m512 gami = _mm512_load_ps(&pgami[0]);
                         const __m512 pi24 = _mm512_set1_ps(2.467401100272339654708622749969f);
                         const __m512 _2pi = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 c0   = _mm512_set1_ps(0.8905f);
                         register __m512 rat1,arg,sarg,arg2,larg2,k0h,t0,rat;
                         register __m512 rcs,cpsi,cgami,cgams,c2gami,c2gams,spsi;
                         register __m512 x0,x1;
                         k0h   = _mm512_mul_ps(k0,h);
                         t0    = _mm512_mul_ps(_2pi,_mm512_mul_ps(h,h));
                         x0    = _mm512_add_ps(k0h,k0h);
                         spsi  = xsinf(psi);
                         arg   = _mm512_mul_ps(x0,spsi);
                         cpsi  = xcosf(psi);
                         arg2  = _mm512_mul_ps(cpsi,_mm512_mul_ps(k0a,c0));
                         larg  = _mm512_fmadd_ps(arg2,arg2,pi24);
                         sarg  = xsinf(arg);
                         cgams = xcosf(gams);
                         rat   = _mm512_div_ps(sarg,arg);
                         cgami = xcosf(gami);
                         x1    = _mm512_mul_ps(rat,rat);
                         c2gams= _mm512_mul_ps(cgams,cgams);
                         c2gami= _mm512_mul_ps(cgami,cgami);
                         x0    = _mm512_mul_ps(t0,_mm512_mul_ps(c2gams,c2gami));
                         rat1  = _mm512_div_ps(x0,larg);
                         rcs   = _mm512_mul_ps(rat1,x1);
                         return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4345_zmm16r4_u(const float * __restrict   ph,
                                              const float * __restrict   pk0,
                                              const float * __restrict   pk0a,
                                              const float * __restrict   pgams,
                                              const float * __restrict   pgami,
                                              const float * __restrict   ppsi) {

                         register __m512 h    = _mm512_loadu_ps(&ph[0]);
                         register __m512 k0   = _mm512_loadu_ps(&pk0[0]); 
                         register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 psi  = _mm512_loadu_ps(&ppsi[0]);
                         register __m512 gams = _mm512_loadu_ps(&pgams[0]);
                         register __m512 gami = _mm512_loadu_ps(&pgami[0]);
                         const __m512 pi24 = _mm512_set1_ps(2.467401100272339654708622749969f);
                         const __m512 _2pi = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 c0   = _mm512_set1_ps(0.8905f);
                         register __m512 rat1,arg,sarg,arg2,larg2,k0h,t0,rat;
                         register __m512 rcs,cpsi,cgami,cgams,c2gami,c2gams,spsi;
                         register __m512 x0,x1;
                         k0h   = _mm512_mul_ps(k0,h);
                         t0    = _mm512_mul_ps(_2pi,_mm512_mul_ps(h,h));
                         x0    = _mm512_add_ps(k0h,k0h);
                         spsi  = xsinf(psi);
                         arg   = _mm512_mul_ps(x0,spsi);
                         cpsi  = xcosf(psi);
                         arg2  = _mm512_mul_ps(cpsi,_mm512_mul_ps(k0a,c0));
                         larg  = _mm512_fmadd_ps(arg2,arg2,pi24);
                         sarg  = xsinf(arg);
                         cgams = xcosf(gams);
                         rat   = _mm512_div_ps(sarg,arg);
                         cgami = xcosf(gami);
                         x1    = _mm512_mul_ps(rat,rat);
                         c2gams= _mm512_mul_ps(cgams,cgams);
                         c2gami= _mm512_mul_ps(cgami,cgami);
                         x0    = _mm512_mul_ps(t0,_mm512_mul_ps(c2gams,c2gami));
                         rat1  = _mm512_div_ps(x0,larg);
                         rcs   = _mm512_mul_ps(rat1,x1);
                         return (rcs);
                }


                  /*
                        Backscattering From a Perfectly Conducting Cylinder With Flat Ends.
                        Helper functions, M1,M2 for the main formula 4.3-48
                        Formula 4.3-50

                   */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 M1_f4350_zmm16r4(const __m512 psi) {

                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,M1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          M1   = _mm512_mul_ps(c2,_mm512_add_ps(x0,inv2));
                          return (M1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 M1_f4350_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,M1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          M1   = _mm512_mul_ps(c2,_mm512_add_ps(x0,inv2));
                          return (M1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 M1_f4350_zmm16r4_u(const float * __restrict  ppsi) {

                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,M1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          M1   = _mm512_mul_ps(c2,_mm512_add_ps(x0,inv2));
                          return (M1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 M2_f4350_zmm16r4(const __m512 psi) {

                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,M2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          M2   = _mm512_mul_ps(c2,_mm512_sub_ps(x0,inv2));
                          return (M2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 M2_f4350_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,M2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          M2   = _mm512_mul_ps(c2,_mm512_sub_ps(x0,inv2));
                          return (M2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 M2_f4350_zmm16r4_u(const float * __restrict  ppsi) {

                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,M2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          M2   = _mm512_mul_ps(c2,_mm512_sub_ps(x0,inv2));
                          return (M2);
                 }


                    /*
                        Backscattering From a Perfectly Conducting Cylinder With Flat Ends.
                        Helper functions, M1,M2 for the main formula 4.3-48
                        Formula 4.3-51

                   */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 N1_f4351_zmm16r4(const __m512 psi) {

                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,N1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          N1   = _mm512_mul_ps(c2,_mm512_sub_ps(n4,_mm512_sub_ps(x0,inv2)));
                          return (N1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 N1_f4351_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,N1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          N1   = _mm512_mul_ps(c2,_mm512_sub_ps(n4,_mm512_sub_ps(x0,inv2)));
                          return (N1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 N1_f4351_zmm16r4_u(const float * __restrict  ppsi) {

                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,N1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          N1   = _mm512_mul_ps(c2,_mm512_sub_ps(n4,_mm512_sub_ps(x0,inv2)));
                          return (N1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 N2_f4351_zmm16r4(const __m512 psi) {

                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,N2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          N2   = _mm512_mul_ps(c2,_mm512_add_ps(n4,_mm512_add_ps(x0,inv2)));
                          return (N2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 N2_f4351_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,N2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          N2   = _mm512_mul_ps(c2,_mm512_add_ps(n4,_mm512_add_ps(x0,inv2)));
                          return (N2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 N2_f4351_zmm16r4_u(const float * __restrict  ppsi) {

                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 n4 = _mm512_set1_ps(-4.0f);
                          const __m512 n1 = _mm512_set1_ps(-1.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          const __m512 c3 = _mm512_set1_ps(0.666666666666666666666666666667f);
                          const __m512 _2 = _mm512_set1_ps(2.0f);
                          register __m512 inv1,inv2,N2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_ps(_2,psi,PI);
                          carg1= _mm512_add_ps(c1,carg1);
                          arg2 = _mm512_mul_ps(c3,x0);
                          inv1 = _mm512_rcp14_ps(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_ps(c1,carg2);
                          inv2 = _mm512_rcp14_ps(x1);
                          x0   = _mm512_mul_ps(n1,inv1);
                          N2   = _mm512_mul_ps(c2,_mm512_add_ps(n4,_mm512_add_ps(x0,inv2)));
                          return (N2);
                }


                   /*
                        Backscattering From a Perfectly Conducting Cylinder With Flat Ends.
                        Helper functions, M1,M2 for the main formula 4.3-48
                        Formula 4.3-52

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G_f4352_zmm16r4(const __m512 psi) {

                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 _2 = _mm512_set1_ps(-2.0f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          register __m512 G,inv,arg,carg,x0;
                          arg = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_ps(c1,carg);
                          inv = _mm512_rcp14_ps(x0);
                          G   = _mm512_mul_ps(c2,_mm512_sub_ps(n2,inv));
                          return (G);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G_f4352_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 _2 = _mm512_set1_ps(-2.0f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          register __m512 G,inv,arg,carg,x0;
                          arg = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_ps(c1,carg);
                          inv = _mm512_rcp14_ps(x0);
                          G   = _mm512_mul_ps(c2,_mm512_sub_ps(n2,inv));
                          return (G);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 G_f4352_zmm16r4_u(const float * __restrict  ppsi) {

                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 _2 = _mm512_set1_ps(-2.0f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          register __m512 G,inv,arg,carg,x0;
                          arg = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_ps(c1,carg);
                          inv = _mm512_rcp14_ps(x0);
                          G   = _mm512_mul_ps(c2,_mm512_sub_ps(n2,inv));
                          return (G);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F_f4352_zmm16r4(const __m512 psi) {

                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 _2 = _mm512_set1_ps(-2.0f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          register __m512 F,inv,arg,carg,x0;
                          arg = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_ps(c1,carg);
                          inv = _mm512_rcp14_ps(x0);
                          G   = _mm512_mul_ps(c2,_mm512_add_ps(n2,inv));
                          return (F);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F_f4352_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512 psi = _mm512_load_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 _2 = _mm512_set1_ps(-2.0f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          register __m512 F,inv,arg,carg,x0;
                          arg = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_ps(c1,carg);
                          inv = _mm512_rcp14_ps(x0);
                          G   = _mm512_mul_ps(c2,_mm512_add_ps(n2,inv));
                          return (F);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 F_f4352_zmm16r4_u(const float * __restrict ppsi) {

                          register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
                          const __m512 c0 = _mm512_set1_ps(0.333333333333333333333333333333333333333333f);
                          const __m512 _2 = _mm512_set1_ps(-2.0f);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          const __m512 c2 = _mm512_set1_ps(0.577350269189625764509148780502f);
                          register __m512 F,inv,arg,carg,x0;
                          arg = _mm512_mul_ps(_mm512_mul_ps(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_ps(c1,carg);
                          inv = _mm512_rcp14_ps(x0);
                          G   = _mm512_mul_ps(c2,_mm512_add_ps(n2,inv));
                          return (F);
                  }


                    /*
                           Scattering From Cylinder Near the Specular Direction.
                           Formula 4.3-53
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4353_zmm16r4(const __m512 k0a,
                                            const __m512 k0,
                                            const __m512 h,
                                            const __m512 phi,
                                            const __m512 psii,
                                            const __m512 psis) {

                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,trm2,trm3;
                          register __m512 cphi,cpsis,c2psis,cpsii,c2psii;
                          register __m512 spsii,spsis,arg,sarg,x0,x1;
                          x0    = _mm512_mul_ps(h,h);
                          x1    = _mm512_mul_ps(c1,phi);
                          trm1  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,x0));
                          cpsii = xcosf(psi);
                          cphi  = xcosf(x1);
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          x0    = _mm512_add_ps(spsii,spsis);
                          c2psis= _mm512_mul_ps(cpsis,cpsis);
                          arg   = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                          x1    = _mm512_mul_ps(c2psis,cphi);
                          sarg  = xsinf(arg);
                          trm2  = _mm512_div_ps(x1,cpsii);
                          trm3  = _mm512_div_ps(sarg,arg);
                          x1    = _mm512_mul_ps(trm1,trm2);
                          x0    = _mm512_mul_ps(trm3,trm3)
                          rcs   = _mm512_mul_ps(x1,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4353_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) ph,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsis) {

                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]);
                          register __m512 h    = _mm512_load_ps(&ph[0]);
                          register __m512 phi  = _mm512_load_ps(&pphi[0]);
                          register __m512 psii = _mm512_load_ps(&ppsii[0]);
                          register __m512 psis = _mm512_load_ps(&ppsis[0]);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,trm2,trm3;
                          register __m512 cphi,cpsis,c2psis,cpsii,c2psii;
                          register __m512 spsii,spsis,arg,sarg,x0,x1;
                          x0    = _mm512_mul_ps(h,h);
                          x1    = _mm512_mul_ps(c1,phi);
                          trm1  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,x0));
                          cpsii = xcosf(psi);
                          cphi  = xcosf(x1);
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          x0    = _mm512_add_ps(spsii,spsis);
                          c2psis= _mm512_mul_ps(cpsis,cpsis);
                          arg   = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                          x1    = _mm512_mul_ps(c2psis,cphi);
                          sarg  = xsinf(arg);
                          trm2  = _mm512_div_ps(x1,cpsii);
                          trm3  = _mm512_div_ps(sarg,arg);
                          x1    = _mm512_mul_ps(trm1,trm2);
                          x0    = _mm512_mul_ps(trm3,trm3)
                          rcs   = _mm512_mul_ps(x1,x0);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4353_zmm16r4_u(const float * __restrict  pk0a,
                                              const float * __restrict  pk0,
                                              const float * __restrict  ph,
                                              const float * __restrict  pphi,
                                              const float * __restrict  ppsii,
                                              const float * __restrict  ppsis) {

                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                          register __m512 h    = _mm512_loadu_ps(&ph[0]);
                          register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                          register __m512 psii = _mm512_loadu_ps(&ppsii[0]);
                          register __m512 psis = _mm512_loadu_ps(&ppsis[0]);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,trm2,trm3;
                          register __m512 cphi,cpsis,c2psis,cpsii,c2psii;
                          register __m512 spsii,spsis,arg,sarg,x0,x1;
                          x0    = _mm512_mul_ps(h,h);
                          x1    = _mm512_mul_ps(c1,phi);
                          trm1  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,x0));
                          cpsii = xcosf(psi);
                          cphi  = xcosf(x1);
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          x0    = _mm512_add_ps(spsii,spsis);
                          c2psis= _mm512_mul_ps(cpsis,cpsis);
                          arg   = _mm512_mul_ps(k0,_mm512_mul_ps(x0,h));
                          x1    = _mm512_mul_ps(c2psis,cphi);
                          sarg  = xsinf(arg);
                          trm2  = _mm512_div_ps(x1,cpsii);
                          trm3  = _mm512_div_ps(sarg,arg);
                          x1    = _mm512_mul_ps(trm1,trm2);
                          x0    = _mm512_mul_ps(trm3,trm3)
                          rcs   = _mm512_mul_ps(x1,x0);
                          return (rcs);
                 }


                   /*

                            Specular direction -- RCS.
                            Formula 4.3-54
                       */

                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4354_zmm16r4(const __m512 k0a,
                                            const __m512 h,
                                            const __m512 psii,
                                            const __m512 phi) {

                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,phi2;
                          register __m512 h2,cpsii,cphi,x0;
                          h2    = _mm512_mul_ps(h,h);
                          phi2  = _mm512_mul_ps(c1,phi);
                          trm1  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          cpsii = xcosf(psii);
                          x0    = _mm512_mul_ps(trm1,cpsii);
                          cphi  = xcosf(phi2);
                          rcs   = _mm512_mul_ps(x0,cphi);
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4354_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ph,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512  k0a  = _mm512_load_ps(&pk0a[0]);
                          register __m512  h    = _mm512_load_ps(&ph[0]);
                          register __m512  psii = _mm512_load_ps(&ppsii[0]);
                          register __m512  phi  = _mm512_load_ps(&pphi[0]);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,phi2;
                          register __m512 h2,cpsii,cphi,x0;
                          h2    = _mm512_mul_ps(h,h);
                          phi2  = _mm512_mul_ps(c1,phi);
                          trm1  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          cpsii = xcosf(psii);
                          x0    = _mm512_mul_ps(trm1,cpsii);
                          cphi  = xcosf(phi2);
                          rcs   = _mm512_mul_ps(x0,cphi);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4354_zmm16r4_u(const float * __restrict  pk0a,
                                              const float * __restrict  ph,
                                              const float * __restrict  ppsii,
                                              const float * __restrict  pphi) {

                          register __m512  k0a  = _mm512_loadu_ps(&pk0a[0]);
                          register __m512  h    = _mm512_loadu_ps(&ph[0]);
                          register __m512  psii = _mm512_loadu_ps(&ppsii[0]);
                          register __m512  phi  = _mm512_loadu_ps(&pphi[0]);
                          const __m512 c1 = _mm512_set1_ps(0.5f);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,phi2;
                          register __m512 h2,cpsii,cphi,x0;
                          h2    = _mm512_mul_ps(h,h);
                          phi2  = _mm512_mul_ps(c1,phi);
                          trm1  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          cpsii = xcosf(psii);
                          x0    = _mm512_mul_ps(trm1,cpsii);
                          cphi  = xcosf(phi2);
                          rcs   = _mm512_mul_ps(x0,cphi);
                          return (rcs);
                 }


                  /*

                         Backscattering direction -- RCS for incidence angles
                         near broadside.
                         Formula 4.3-54
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4354_zmm16r4(const __m512 k0a,
                                            const __m512 h,
                                            const __m512 k0,
                                            const __m512 psii) {

                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,trm2,cpsii,spsii;
                          register __m512 x0,x1,k0h,h2,arg,sarg;
                          k0h  = _mm512_mul_ps(k0,h);
                          h2   = _mm512_mul_ps(h,h);
                          x0   = _mm512_add_ps(k0h,k0h); 
                          x1   = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          cpsii= xcosf(psi);
                          spsii= xsinf(psi);
                          trm1 = _mm512_mul_ps(x1,cpsii);
                          arg  = _mm512_mul_ps(x0,spsii);
                          sarg = xsinf(arg);
                          x0   = _mm512_div_ps(sarg,arg);
                          trm2 = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4354_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ph,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) ppsii) {

                          register __m512  k0a  = _mm512_load_ps(&pk0a[0]);
                          register __m512  h    = _mm512_load_ps(&ph[0]);
                          register __m512  psii = _mm512_load_ps(&ppsii[0]);
                          register __m512  k0   = _mm512_load_ps(&pk0[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,trm2,cpsii,spsii;
                          register __m512 x0,x1,k0h,h2,arg,sarg;
                          k0h  = _mm512_mul_ps(k0,h);
                          h2   = _mm512_mul_ps(h,h);
                          x0   = _mm512_add_ps(k0h,k0h); 
                          x1   = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          cpsii= xcosf(psi);
                          spsii= xsinf(psi);
                          trm1 = _mm512_mul_ps(x1,cpsii);
                          arg  = _mm512_mul_ps(x0,spsii);
                          sarg = xsinf(arg);
                          x0   = _mm512_div_ps(sarg,arg);
                          trm2 = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4354_zmm16r4_u(const float * __restrict  pk0a,
                                              const float * __restrict  ph,
                                              const float * __restrict  pk0,
                                              const float * __restrict  ppsii) {

                          register __m512  k0a  = _mm512_loadu_ps(&pk0a[0]);
                          register __m512  h    = _mm512_loadu_ps(&ph[0]);
                          register __m512  psii = _mm512_loadu_ps(&ppsii[0]);
                          register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,trm1,trm2,cpsii,spsii;
                          register __m512 x0,x1,k0h,h2,arg,sarg;
                          k0h  = _mm512_mul_ps(k0,h);
                          h2   = _mm512_mul_ps(h,h);
                          x0   = _mm512_add_ps(k0h,k0h); 
                          x1   = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          cpsii= xcosf(psi);
                          spsii= xsinf(psi);
                          trm1 = _mm512_mul_ps(x1,cpsii);
                          arg  = _mm512_mul_ps(x0,spsii);
                          sarg = xsinf(arg);
                          x0   = _mm512_div_ps(sarg,arg);
                          trm2 = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(trm1,trm2);
                          return (rcs);
                }


                 /*

                        Broadside (psi == 0) RCS.
                        Formula 4.3-56
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4356_zmm16r4(const __m512 k0a,
                                            const __m512 h) {

                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,h2;
                          h2 = _mm512_mul_ps(h,h);
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          return (rcs); 
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4356_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(64) ph) {

                          register __m512  k0a  = _mm512_load_ps(&pk0a[0]);
                          register __m512  h    = _mm512_load_ps(&ph[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,h2;
                          h2 = _mm512_mul_ps(h,h);
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          return (rcs); 
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4356_zmm16r4_u(const float * __restrict  pk0a,
                                            const float * __restrict  ph) {

                          register __m512  k0a  = _mm512_loadu_ps(&pk0a[0]);
                          register __m512  h    = _mm512_loadu_ps(&ph[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,h2;
                          h2 = _mm512_mul_ps(h,h);
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,h2));
                          return (rcs); 
                }


                  /*
                       Elliptical cylinders.
                   */


                   /*
                         Low-frequency approximations (k0a<0.5, k0b<0.5)
                         TM-case,formula 4.4-11
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4411_zmm16r4(const __m512 a,
                                         const __m512 b,
                                         const __m512 k0,
                                         __m512 * __restrict TMr,
                                         __m512 * __restrict TMi) {

                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         const __m512 imn  = _mm512_set1_ps(-1.57079632679489661923132169164f);
                         const __m512 imp  = _mm512_set1_ps(1.57079632679489661923132169164f);
                         const __m512 c0   = _mm512_set1_ps(0.8905f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 ab2,c0k0,arg,larg;
                         register __m512 invr,invi;
                         ab2  = _mm512_mul_ps(_mm512_add_ps(a,b),hlf);
                         c0k0 = _mm512_mul_ps(c0,k0);
                         arg  = _mm512_mul_ps(ab2,c0k0);
                         larg = xlogf(arg);
                         cdiv_zmm16r4(_1,_1,larg,imn,&invr,&invi);
                         *TMr = _mm512_mul_ps(imp,invr);
                         *TMi = _mm512_mul_ps(imp,invi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4411_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pb,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         float * __restrict __ATTR_ALIGN__(64) TMr,
                                         float * __restrict __ATTR_ALIGN__(64) TMi) {

                         register __m512 a = _mm512_load_ps(&pa[0]);
                         register __m512 b = _mm512_load_ps(&pb[0]);
                         register __m512 k0= _mm512_load_ps(&pk0[0]);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         const __m512 imn  = _mm512_set1_ps(-1.57079632679489661923132169164f);
                         const __m512 imp  = _mm512_set1_ps(1.57079632679489661923132169164f);
                         const __m512 c0   = _mm512_set1_ps(0.8905f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 ab2,c0k0,arg,larg;
                         register __m512 invr,invi;
                         ab2  = _mm512_mul_ps(_mm512_add_ps(a,b),hlf);
                         c0k0 = _mm512_mul_ps(c0,k0);
                         arg  = _mm512_mul_ps(ab2,c0k0);
                         larg = xlogf(arg);
                         cdiv_zmm16r4(_1,_1,larg,imn,&invr,&invi);
                         _mm512_store_ps(&TMr[0], _mm512_mul_ps(imp,invr));
                         _mm512_store_ps(&TMi[0], _mm512_mul_ps(imp,invi));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4411_zmm16r4_u(const float * __restrict  pa,
                                           const float * __restrict  pb,
                                           const float * __restrict  pk0,
                                           float * __restrict  TMr,
                                           float * __restrict  TMi) {

                         register __m512 a = _mm512_loadu_ps(&pa[0]);
                         register __m512 b = _mm512_loadu_ps(&pb[0]);
                         register __m512 k0= _mm512_loadu_ps(&pk0[0]);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         const __m512 imn  = _mm512_set1_ps(-1.57079632679489661923132169164f);
                         const __m512 imp  = _mm512_set1_ps(1.57079632679489661923132169164f);
                         const __m512 c0   = _mm512_set1_ps(0.8905f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
                         register __m512 ab2,c0k0,arg,larg;
                         register __m512 invr,invi;
                         ab2  = _mm512_mul_ps(_mm512_add_ps(a,b),hlf);
                         c0k0 = _mm512_mul_ps(c0,k0);
                         arg  = _mm512_mul_ps(ab2,c0k0);
                         larg = xlogf(arg);
                         cdiv_zmm16r4(_1,_1,larg,imn,&invr,&invi);
                         _mm512_storeu_ps(&TMr[0], _mm512_mul_ps(imp,invr));
                         _mm512_storeu_ps(&TMi[0], _mm512_mul_ps(imp,invi));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4412_zmm16r4(const __m512 k0a,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 phi1,
                                         const __m512 phi2,
                                         __m512 * __restrict TEr,
                                         __m512 * __restrict TEi) {

                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 k0a2,ba,cphi1,sphi1,trm1,trm2,_1ba,x0,x1;
                        register __m512 cphi2,sphi2;
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        ba    = _mm512_div_ps(b,a);
                        cphi1 = xcosf(phi1);
                        _1ba  = _mm512_add_ps(_1,ba);
                        sphi1 = xsinf(phi1);
                        x0    = _mm512_mul_ps(pi4,k0a2);
                        cphi2 = xcosf(phi2);
                        x1    = _mm512_add_ps(ba,_1ba);
                        sphi2 = xsinf(phi2);
                        trm1  = _mm512_mul_ps(x0,x1);
                        x0    = _mm512_fmadd_ps(cphi2,cphi1,_mm512_mul_ps(sphi2,sphi1));
                        trm2  = _mm512_mul_ps(ba,x0);
                        x1    = _mm512_mul_ps(trm1,trm2);
                        *TEr  = nIi;
                        *TEi  = _mm512_mul_ps(nIi,x1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4412_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pb,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           float * __restrict __ATTR_ALIGN__(64) TEr,
                                           float * __restrict __ATTR_ALIGN__(64) TEi) {

                        register __m512 a    = _mm512_load_ps(&pa[0]);
                        register __m512 b    = _mm512_load_ps(&pb[0]);
                        register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        register __m512 phi1 = _mm512_load_ps(&phi1[0]);
                        register __m512 phi2 = _mm512_load_ps(&phi2[0]);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 k0a2,ba,cphi1,sphi1,trm1,trm2,_1ba,x0,x1;
                        register __m512 cphi2,sphi2;
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        ba    = _mm512_div_ps(b,a);
                        cphi1 = xcosf(phi1);
                        _1ba  = _mm512_add_ps(_1,ba);
                        sphi1 = xsinf(phi1);
                        x0    = _mm512_mul_ps(pi4,k0a2);
                        cphi2 = xcosf(phi2);
                        x1    = _mm512_add_ps(ba,_1ba);
                        sphi2 = xsinf(phi2);
                        trm1  = _mm512_mul_ps(x0,x1);
                        x0    = _mm512_fmadd_ps(cphi2,cphi1,_mm512_mul_ps(sphi2,sphi1));
                        trm2  = _mm512_mul_ps(ba,x0);
                        x1    = _mm512_mul_ps(trm1,trm2);
                        _mm512_store_ps(&TEr[0] ,nIi);
                        _mm512_store_ps(&TEi[0] ,_mm512_mul_ps(nIi,x1));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4412_zmm16r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  pa,
                                           const float * __restrict  pb,
                                           const float * __restrict  pphi1,
                                           const float * __restrict  pphi2,
                                           __m512 * __restrict TEr,
                                           __m512 * __restrict TEi) {

                        register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        register __m512 b    = _mm512_loadu_ps(&pb[0]);
                        register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512 phi1 = _mm512_loadu_ps(&phi1[0]);
                        register __m512 phi2 = _mm512_loadu_ps(&phi2[0]);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        register __m512 k0a2,ba,cphi1,sphi1,trm1,trm2,_1ba,x0,x1;
                        register __m512 cphi2,sphi2;
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        ba    = _mm512_div_ps(b,a);
                        cphi1 = xcosf(phi1);
                        _1ba  = _mm512_add_ps(_1,ba);
                        sphi1 = xsinf(phi1);
                        x0    = _mm512_mul_ps(pi4,k0a2);
                        cphi2 = xcosf(phi2);
                        x1    = _mm512_add_ps(ba,_1ba);
                        sphi2 = xsinf(phi2);
                        trm1  = _mm512_mul_ps(x0,x1);
                        x0    = _mm512_fmadd_ps(cphi2,cphi1,_mm512_mul_ps(sphi2,sphi1));
                        trm2  = _mm512_mul_ps(ba,x0);
                        x1    = _mm512_mul_ps(trm1,trm2);
                        _mm512_store_ps(&TEr[0], nIi);
                        _mm512_store_ps(&TEi[0], _mm512_mul_ps(nIi,x1));
                }


                 /*
                       TM-case, RCS.
                       Formula 4.4-13
                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4413_zmm16r4(const __m512 a,
                                            const __m512 b,
                                            const __m512 k0) {

                          const __m512 pi2 = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(0.8905f);
                          const __m512 pi24= _mm512_set1_ps(2.467401100272339654708622749969f);
                          register __m512 rcs,abh,k0abh,num,sqr1,sqr2,c0k0,arg,larg,x0,den,x1;
                          abh = _mm512_mul_ps(_mm512_add_ps(a,b),hlf);
                          c0k0= _mm512_mul_ps(c0,k0);
                          num = _mm512_mul_ps(pi2,abh);
                          arg = _mm512_mul_ps(c0k0,abh);
                          larg= xlogf(arg);
                          x0  = _mm512_fmadd_ps(larg,larg,pi24);
                          sqr1= _mm512_sqrt_ps(_mm512_mul_ps(k0,abh));
                          sqr2= _mm512_sqrt_ps(x0);
                          den = _mm512_mul_ps(sqr1,sqr2);
                          x1  = _mm512_mul_ps(den,den);
                          rcs = _mm512_div_ps(num,x1);
                          return (rcs);
                }
                                            


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4413_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pa,
                                              const float * __restrict __ATTR_ALIGN__(64)  pb,
                                              const float * __restrict __ATTR_ALIGN__(64)  pk0) {

                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]);
                          const __m512 pi2 = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(0.8905f);
                          const __m512 pi24= _mm512_set1_ps(2.467401100272339654708622749969f);
                          register __m512 rcs,abh,k0abh,num,sqr1,sqr2,c0k0,arg,larg,x0,den,x1;
                          abh = _mm512_mul_ps(_mm512_add_ps(a,b),hlf);
                          c0k0= _mm512_mul_ps(c0,k0);
                          num = _mm512_mul_ps(pi2,abh);
                          arg = _mm512_mul_ps(c0k0,abh);
                          larg= xlogf(arg);
                          x0  = _mm512_fmadd_ps(larg,larg,pi24);
                          sqr1= _mm512_sqrt_ps(_mm512_mul_ps(k0,abh));
                          sqr2= _mm512_sqrt_ps(x0);
                          den = _mm512_mul_ps(sqr1,sqr2);
                          x1  = _mm512_mul_ps(den,den);
                          rcs = _mm512_div_ps(num,x1);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4413_zmm16r4_u(const float * __restrict   pa,
                                              const float * __restrict   pb,
                                              const float * __restrict   pk0) {

                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                          const __m512 pi2 = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(0.8905f);
                          const __m512 pi24= _mm512_set1_ps(2.467401100272339654708622749969f);
                          register __m512 rcs,abh,k0abh,num,sqr1,sqr2,c0k0,arg,larg,x0,den,x1;
                          abh = _mm512_mul_ps(_mm512_add_ps(a,b),hlf);
                          c0k0= _mm512_mul_ps(c0,k0);
                          num = _mm512_mul_ps(pi2,abh);
                          arg = _mm512_mul_ps(c0k0,abh);
                          larg= xlogf(arg);
                          x0  = _mm512_fmadd_ps(larg,larg,pi24);
                          sqr1= _mm512_sqrt_ps(_mm512_mul_ps(k0,abh));
                          sqr2= _mm512_sqrt_ps(x0);
                          den = _mm512_mul_ps(sqr1,sqr2);
                          x1  = _mm512_mul_ps(den,den);
                          rcs = _mm512_div_ps(num,x1);
                          return (rcs);
                }


                    /*
                         High frequency approximations (k0a>5, k0b>5)
                         TM-case, formula 4.4-15
                      */

#include "GMS_simd_utils.hpp"


                    /*
                        Helper function for testing the condition of high-frequency limit.
                        Page. 322.

                     */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __mmask16 
                   TM_f4415_helper_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 phi1,
                                           const __m512 phi2,
                                           const __m512 b) {
                      const __m512 c0 = _mm512_set1_ps(0.166666666666666666666666666667f);
                      register __m512 a2,b2,sphi1,cphi1,trm1,trm2,root6;
                      register __m512 k02,absp,sphi1s,cphi1s,k0a2,k0b2,x0;
                      __mmask16 m;
                      k02  = _mm512_mul_ps(k0,k0);
                      a2   = _mm512_mul_ps(a,a);
                      k0a2 = _mm512_mul_ps(k02,a2);
                      b2   = _mm512_mul_ps(b,b);
                      k0b2 = _mm512_mul_ps(k02,b2);
                      cphi1= xcosf(phi1);
                      absp = _mm512_abs_ps(_mm512_sub_ps(phi2,phi1));
                      cphi1s = _mm512_mul_ps(cphi1,cphi1)
                      sphi1= xsinf(phi1);
                      trm1 = _mm512_sub_ps(PI,absp);
                      sphi1s = _mm512_mul_ps(sphi1,sphi1);
                      trm2 = _mm512_fmadd_ps(k02a2,sphi1s,_mm512_mul_ps(k02b2,cphi1s));
                      x0   = _mm512_pow_ps(trm2,c0);
                      root6= _mm512_rcp14_ps(x0);
                      m    = _mm512_cmp_mask_ps(trm1,root6,_CMP_GT_OQ);
                      return (m);
                }

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4415_zmm16r4(const __m512 phi1,
                                         const __m512 phi2,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 k0,
                                         __m512 * __restrict TMr,
                                         __m512 * __restrict TMi,
                                         bool & status) {

                        using namespace gms::math;
                        __mmask16 m = TM_f4415_helper_zmm16r4(k0,a,phi1,phi2,b);
                        if(!m) {
                           status = false;
                           return;
                        }
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 ip4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 arg1,arg2,carg1,carg2,sarg2,sqr1,ear,eai,cer,cei,trm1;
                        register __m512 f,rho,a2b2,a2,b2,b2a2,k0a,cphi2,cphi1,sphi1,sphi2,frat;
                        register __m512 cphis,sphis,rhod,rhorat,x0,x1,tmp1,tmp2,b2a2s,carg2s,sarg2s;
                        arg1  = _mm512_mul_ps(_mm512_sub_ps(phi2,phi1),hlf);
                        a2    = _mm512_mul_ps(a,a);
                        b2    = _mm512_mul_ps(b,b);
                        k0a   = _mm512_mul_ps(k0,a);
                        arg2  = _mm512_mul_ps(_mm512_add_ps(phi2,phi1),hlf);
                        carg1 = xcosf(arg1);
                        a2b2  = _mm512_mul_ps(a2,b2);
                        b2a2  = _mm512_div_ps(b2,a2);
                        cphi1 = xcosf(phi1);
                        sphi1 = xsinf(phi1);
                        trm1  = _mm512_sqrt_ps(_mm512_mul_ps(PI,carg1));
                        cphi2 = xcosf(phi2);
                        sphi2 = xsinf(phi2);
                        cphis = _mm512_add_ps(cphi1,cphi2);
                        carg2 = xcosf(arg2);
                        sphis = _mm512_add_ps(sphi1,sphi2);
                        sarg2 = xsinf(arg2);
                        x0    = _mm512_mul_ps(carg2,carg2);
                        x1    = _mm512_mul_ps(sarg2,sarg2);
                        rhod  = _mm512_fmadd_ps(a2,x0,_mm512_mul_ps(b2,x1));
                        b2a2s = _mm512_mul_ps(b2a2,sphis);
                        tmp1  = _mm512_pow_ps(rhod,_mm512_set1_ps(1.5f));
                        rhorat= _mm512_div_ps(a2b2,tmp1);
                        x0    = _mm512_fmadd_ps(sarg2,b2a2s,carg2);
                        carg2s= _mm512_mul_ps(carg2,carg2)
                        tmp2  = _mm512_mul_ps(cphis,x0);
                        sarg2s= _mm512_mul_ps(sarg2,sarg2);
                        x1    = _mm512_fmadd_ps(b2a2,sarg2s,carg2s);
                        tmp1  = _mm512_sqrt_ps(x1);
                        frat  = _mm512_div_ps(tmp2,tmp1);
                        trm1  = zmm16r4_negate(trm1);
                        ear   = _mm512_add_ps(nIr,ip4);
                        x0    = _mm512_mul_ps(_mm512_sqrt_ps(_mm512_mul_ps(k0,rhorat)),hlf);
                        eai   = _mm512_mul_ps(nIi,_mm512_mul_ps(k0a,frat));
                        eai   = _mm512_add_ps(eai,ip4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        x1    = _mm512_mul_ps(trm1,x0);
                        *TMr = _mm512_mul_ps(x1,cer);
                        *TMi = _mm512_mul_ps(x1,cei);
                        status = true;
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4415_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pb,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         float * __restrict __ATTR_ALIGN__(64) TMr,
                                         float * __restrict __ATTR_ALIGN__(64) TMi,
                                         bool & status) {

                        using namespace gms::math;

                        register __m512 phi1 = _mm512_load_ps(&phi1[0]);
                        register __m512 phi2 = _mm512_load_ps(&phi2[0]);
                        register __m512 a    = _mm512_load_ps(&pa[0]);
                        register __m512 b    = _mm512_load_ps(&pb[0]);
                        register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        __mmask16 m = TM_f4415_helper_zmm16r4(k0,a,phi1,phi2,b);
                        if(!m) {
                           status = false;
                           return;
                        }
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 ip4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 arg1,arg2,carg1,carg2,sarg2,sqr1,ear,eai,cer,cei,trm1;
                        register __m512 f,rho,a2b2,a2,b2,b2a2,k0a,cphi2,cphi1,sphi1,sphi2,frat;
                        register __m512 cphis,sphis,rhod,rhorat,x0,x1,tmp1,tmp2,b2a2s,carg2s,sarg2s;
                        arg1  = _mm512_mul_ps(_mm512_sub_ps(phi2,phi1),hlf);
                        a2    = _mm512_mul_ps(a,a);
                        b2    = _mm512_mul_ps(b,b);
                        k0a   = _mm512_mul_ps(k0,a);
                        arg2  = _mm512_mul_ps(_mm512_add_ps(phi2,phi1),hlf);
                        carg1 = xcosf(arg1);
                        a2b2  = _mm512_mul_ps(a2,b2);
                        b2a2  = _mm512_div_ps(b2,a2);
                        cphi1 = xcosf(phi1);
                        sphi1 = xsinf(phi1);
                        trm1  = _mm512_sqrt_ps(_mm512_mul_ps(PI,carg1));
                        cphi2 = xcosf(phi2);
                        sphi2 = xsinf(phi2);
                        cphis = _mm512_add_ps(cphi1,cphi2);
                        carg2 = xcosf(arg2);
                        sphis = _mm512_add_ps(sphi1,sphi2);
                        sarg2 = xsinf(arg2);
                        x0    = _mm512_mul_ps(carg2,carg2);
                        x1    = _mm512_mul_ps(sarg2,sarg2);
                        rhod  = _mm512_fmadd_ps(a2,x0,_mm512_mul_ps(b2,x1));
                        b2a2s = _mm512_mul_ps(b2a2,sphis);
                        tmp1  = _mm512_pow_ps(rhod,_mm512_set1_ps(1.5f));
                        rhorat= _mm512_div_ps(a2b2,tmp1);
                        x0    = _mm512_fmadd_ps(sarg2,b2a2s,carg2);
                        carg2s= _mm512_mul_ps(carg2,carg2)
                        tmp2  = _mm512_mul_ps(cphis,x0);
                        sarg2s= _mm512_mul_ps(sarg2,sarg2);
                        x1    = _mm512_fmadd_ps(b2a2,sarg2s,carg2s);
                        tmp1  = _mm512_sqrt_ps(x1);
                        frat  = _mm512_div_ps(tmp2,tmp1);
                        trm1  = zmm16r4_negate(trm1);
                        ear   = _mm512_add_ps(nIr,ip4);
                        x0    = _mm512_mul_ps(_mm512_sqrt_ps(_mm512_mul_ps(k0,rhorat)),hlf);
                        eai   = _mm512_mul_ps(nIi,_mm512_mul_ps(k0a,frat));
                        eai   = _mm512_add_ps(eai,ip4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        x1    = _mm512_mul_ps(trm1,x0);
                        _mm512_store_ps(&TMr[0] ,_mm512_mul_ps(x1,cer));
                        _mm512_store_ps(&TMi[0] ,_mm512_mul_ps(x1,cei));
                        status = true;
                 }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4415_zmm16r4_u(const float * __restrict  pphi1,
                                         const float * __restrict  pphi2,
                                         const float * __restrict  pa,
                                         const float * __restrict  pb,
                                         const float * __restrict  pk0,
                                         float * __restrict  TMr,
                                         float * __restrict  TMi,
                                         bool & status) {

                        using namespace gms::math;
                        register __m512 phi1 = _mm512_loadu_ps(&phi1[0]);
                        register __m512 phi2 = _mm512_loadu_ps(&phi2[0]);
                        register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        register __m512 b    = _mm512_loadu_ps(&pb[0]);
                        register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        __mmask16 m = TM_f4415_helper_zmm16r4(k0,a,phi1,phi2,b);
                        if(!m) {
                           status = false;
                           return;
                        }
                        const __m512 hlf = _mm512_set1_ps(0.5f);
                        const __m512 ip4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 arg1,arg2,carg1,carg2,sarg2,sqr1,ear,eai,cer,cei,trm1;
                        register __m512 f,rho,a2b2,a2,b2,b2a2,k0a,cphi2,cphi1,sphi1,sphi2,frat;
                        register __m512 cphis,sphis,rhod,rhorat,x0,x1,tmp1,tmp2,b2a2s,carg2s,sarg2s;
                        arg1  = _mm512_mul_ps(_mm512_sub_ps(phi2,phi1),hlf);
                        a2    = _mm512_mul_ps(a,a);
                        b2    = _mm512_mul_ps(b,b);
                        k0a   = _mm512_mul_ps(k0,a);
                        arg2  = _mm512_mul_ps(_mm512_add_ps(phi2,phi1),hlf);
                        carg1 = xcosf(arg1);
                        a2b2  = _mm512_mul_ps(a2,b2);
                        b2a2  = _mm512_div_ps(b2,a2);
                        cphi1 = xcosf(phi1);
                        sphi1 = xsinf(phi1);
                        trm1  = _mm512_sqrt_ps(_mm512_mul_ps(PI,carg1));
                        cphi2 = xcosf(phi2);
                        sphi2 = xsinf(phi2);
                        cphis = _mm512_add_ps(cphi1,cphi2);
                        carg2 = xcosf(arg2);
                        sphis = _mm512_add_ps(sphi1,sphi2);
                        sarg2 = xsinf(arg2);
                        x0    = _mm512_mul_ps(carg2,carg2);
                        x1    = _mm512_mul_ps(sarg2,sarg2);
                        rhod  = _mm512_fmadd_ps(a2,x0,_mm512_mul_ps(b2,x1));
                        b2a2s = _mm512_mul_ps(b2a2,sphis);
                        tmp1  = _mm512_pow_ps(rhod,_mm512_set1_ps(1.5f));
                        rhorat= _mm512_div_ps(a2b2,tmp1);
                        x0    = _mm512_fmadd_ps(sarg2,b2a2s,carg2);
                        carg2s= _mm512_mul_ps(carg2,carg2)
                        tmp2  = _mm512_mul_ps(cphis,x0);
                        sarg2s= _mm512_mul_ps(sarg2,sarg2);
                        x1    = _mm512_fmadd_ps(b2a2,sarg2s,carg2s);
                        tmp1  = _mm512_sqrt_ps(x1);
                        frat  = _mm512_div_ps(tmp2,tmp1);
                        trm1  = zmm16r4_negate(trm1);
                        ear   = _mm512_add_ps(nIr,ip4);
                        x0    = _mm512_mul_ps(_mm512_sqrt_ps(_mm512_mul_ps(k0,rhorat)),hlf);
                        eai   = _mm512_mul_ps(nIi,_mm512_mul_ps(k0a,frat));
                        eai   = _mm512_add_ps(eai,ip4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        x1    = _mm512_mul_ps(trm1,x0);
                        _mm512_storeu_ps(&TMr[0] ,_mm512_mul_ps(x1,cer));
                        _mm512_storeu_ps(&TMi[0] ,_mm512_mul_ps(x1,cei));
                        status = true;
                 }


                    /*
                         High frequency approximations (k0a>5, k0b>5)
                         TE-case, formula 4.4-16
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4416_zmm16r4(const __m512 phi1,
                                         const __m512 phi2,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 k0,
                                         __m512 * __restrict TEr,
                                         __m512 * __restrict TEi,
                                         bool & status) {

                        __m512 resr,resi;
                        TM_f4415_zmm16r4(phi1,phi2,a,b,k0,&resr,&resi,status);
                        if(!status) {
                           return;
                        }
                         else {
                           *TEr = _mm512_mul_ps(nIi,resr);
                           *TEi = _mm512_mul_ps(nIi,resi);
                        }
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4416_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pb,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         float * __restrict __ATTR_ALIGN__(64) TEr,
                                         float * __restrict __ATTR_ALIGN__(64) TEi,
                                         bool & status) {

                        
                        TM_f4415_zmm16r4_a(phi1,phi2,a,b,k0,TEr,TEi,status);
                        if(!status) {
                           return;
                        }
                         else {
                           _mm512_store_ps(&TEr[0] ,_mm512_mul_ps(nIi,_mm512_load_ps(&TEr[0])));
                           _mm512_store_ps(&TEi[0] ,_mm512_mul_ps(nIi,_mm512_load_ps(&TEi[0])));
                        }
               }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4416_zmm16r4_u(const float * __restrict  pphi1,
                                         const float * __restrict  pphi2,
                                         const float * __restrict  pa,
                                         const float * __restrict  pb,
                                         const float * __restrict  pk0,
                                         float * __restrict  TEr,
                                         float * __restrict  TEi,
                                         bool & status) {

                       
                        TM_f4415_zmm16r4_u(phi1,phi2,a,b,k0,TEr,TEi,status);
                        if(!status) {
                           return;
                        }
                         else {
                           _mm512_storeu_ps(&TEr[0] ,_mm512_mul_ps(nIi,_mm512_loadu_ps(&TEr[0])));
                           _mm512_storeu_ps(&TEi[0] ,_mm512_mul_ps(nIi,_mm512_loadu_ps(&TEr[0])));
                        }
               }


                 /*

                        Bistatic scattering width.
                        Formula 4.4-19
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4419_zmm16r4(const __m512 phi1,
                                            const __m512 phi2,
                                            const __m512 a,
                                            const __m512 b) {

                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(1.5f);
                          register __m512 rcs,a2,b2,a2b2,num;
                          register __m512 arg,carg,carg2,sarg,sarg2;
                          register __m512 pow32,x0;
                          a2   = _mm512_mul_ps(a,a);
                          arg  = _mm512_mul_ps(_mm512_add_ps(phi2,phi1),hlf);
                          b2   = _mm512_mul_ps(b,b);
                          carg = xcosf(arg);
                          num  = _mm512_mul_ps(PI,_mm512_mul_ps(a2,b2));
                          sarg = xsinf(arg);
                          carg2= _mm512_mul_ps(carg,carg);
                          sarg2= _mm512_mul_ps(sarg,sarg);
                          x0   = _mm512_fmadd_ps(a2,carg2,_mm512_mul_ps(b2,sarg2));
                          pow32= _mm512_pow_ps(x0,c0);
                          rcs  = _mm512_div_ps(num,pow32);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4419_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb) {

                          register __m512 phi1 = _mm512_load_ps(&pphi1[0]);
                          register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(1.5f);
                          register __m512 rcs,a2,b2,a2b2,num;
                          register __m512 arg,carg,carg2,sarg,sarg2;
                          register __m512 pow32,x0;
                          a2   = _mm512_mul_ps(a,a);
                          arg  = _mm512_mul_ps(_mm512_add_ps(phi2,phi1),hlf);
                          b2   = _mm512_mul_ps(b,b);
                          carg = xcosf(arg);
                          num  = _mm512_mul_ps(PI,_mm512_mul_ps(a2,b2));
                          sarg = xsinf(arg);
                          carg2= _mm512_mul_ps(carg,carg);
                          sarg2= _mm512_mul_ps(sarg,sarg);
                          x0   = _mm512_fmadd_ps(a2,carg2,_mm512_mul_ps(b2,sarg2));
                          pow32= _mm512_pow_ps(x0,c0);
                          rcs  = _mm512_div_ps(num,pow32);
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4419_zmm16r4_u(const float * __restrict  pphi1,
                                              const float * __restrict  pphi2,
                                              const float * __restrict  pa,
                                              const float * __restrict  pb) {

                          register __m512 phi1 = _mm512_loadu_ps(&pphi1[0]);
                          register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          const __m512 hlf = _mm512_set1_ps(0.5f);
                          const __m512 c0  = _mm512_set1_ps(1.5f);
                          register __m512 rcs,a2,b2,a2b2,num;
                          register __m512 arg,carg,carg2,sarg,sarg2;
                          register __m512 pow32,x0;
                          a2   = _mm512_mul_ps(a,a);
                          arg  = _mm512_mul_ps(_mm512_add_ps(phi2,phi1),hlf);
                          b2   = _mm512_mul_ps(b,b);
                          carg = xcosf(arg);
                          num  = _mm512_mul_ps(PI,_mm512_mul_ps(a2,b2));
                          sarg = xsinf(arg);
                          carg2= _mm512_mul_ps(carg,carg);
                          sarg2= _mm512_mul_ps(sarg,sarg);
                          x0   = _mm512_fmadd_ps(a2,carg2,_mm512_mul_ps(b2,sarg2));
                          pow32= _mm512_pow_ps(x0,c0);
                          rcs  = _mm512_div_ps(num,pow32);
                          return (rcs);
                 }


                   /*

                          Backscattering width, for phi2 == phi1.
                          Formula 4.4-20
                      */


                     
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4420_zmm16r4(const __m512 a,
                                            const __m512 b,
                                            const __m512 phi) {

                          const __m512 c0  = _mm512_set1_ps(1.5f);
                          register __m512 rcs,a2,b2,a2b2,num;
                          register __m512 carg,carg2,sarg,sarg2;
                          register __m512 pow32,x0;
                          a2   = _mm512_mul_ps(a,a);
                          carg = xcosf(phi);
                          b2   = _mm512_mul_ps(b,b);
                          carg2= _mm512_mul_ps(carg,carg);
                          num  = _mm512_mul_ps(PI,_mm512_mul_ps(a2,b2));
                          sarg = xsinf(phi);
                          sarg2= _mm512_mul_ps(sarg,sarg);
                          x0   = _mm512_fmadd_ps(a2,carg2,_mm512_mul_ps(b2,sarg2));
                          pow32= _mm512_pow_ps(x0,c0);
                          rcs  = _mm512_div_ps(num,pow32);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4420_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512 phi2 = _mm512_load_ps(&pphi[0]);
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          const __m512 c0  = _mm512_set1_ps(1.5f);
                          register __m512 rcs,a2,b2,a2b2,num;
                          register __m512 carg,carg2,sarg,sarg2;
                          register __m512 pow32,x0;
                          a2   = _mm512_mul_ps(a,a);
                          carg = xcosf(phi);
                          b2   = _mm512_mul_ps(b,b);
                          carg2= _mm512_mul_ps(carg,carg);
                          num  = _mm512_mul_ps(PI,_mm512_mul_ps(a2,b2));
                          sarg = xsinf(phi);
                          sarg2= _mm512_mul_ps(sarg,sarg);
                          x0   = _mm512_fmadd_ps(a2,carg2,_mm512_mul_ps(b2,sarg2));
                          pow32= _mm512_pow_ps(x0,c0);
                          rcs  = _mm512_div_ps(num,pow32);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4420_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pphi) {

                          register __m512 phi2 = _mm512_loadu_ps(&pphi[0]);
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          const __m512 c0  = _mm512_set1_ps(1.5f);
                          register __m512 rcs,a2,b2,a2b2,num;
                          register __m512 carg,carg2,sarg,sarg2;
                          register __m512 pow32,x0;
                          a2   = _mm512_mul_ps(a,a);
                          carg = xcosf(phi);
                          b2   = _mm512_mul_ps(b,b);
                          carg2= _mm512_mul_ps(carg,carg);
                          num  = _mm512_mul_ps(PI,_mm512_mul_ps(a2,b2));
                          sarg = xsinf(phi);
                          sarg2= _mm512_mul_ps(sarg,sarg);
                          x0   = _mm512_fmadd_ps(a2,carg2,_mm512_mul_ps(b2,sarg2));
                          pow32= _mm512_pow_ps(x0,c0);
                          rcs  = _mm512_div_ps(num,pow32);
                          return (rcs);
                 }


                   /*
                        Forward scattering pattern and width.
                        Formula 4.4-23 a scattering amplitude

                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __mmask16 
                   T_f4423_helper_zmm16r4( const __m512 k0,
                                           const __m512 a,
                                           const __m512 phi1,
                                           const __m512 phi2,
                                           const __m512 b) {
                      const __m512 c0 = _mm512_set1_ps(0.166666666666666666666666666667f);
                      register __m512 a2,b2,sphi1,cphi1,trm1,trm2,root6;
                      register __m512 k02,absp,sphi1s,cphi1s,k0a2,k0b2,x0;
                      __mmask16 m;
                      k02  = _mm512_mul_ps(k0,k0);
                      a2   = _mm512_mul_ps(a,a);
                      k0a2 = _mm512_mul_ps(k02,a2);
                      b2   = _mm512_mul_ps(b,b);
                      k0b2 = _mm512_mul_ps(k02,b2);
                      cphi1= xcosf(phi1);
                      trm1 = _mm512_sub_ps(phi1,phi2);
                      cphi1s = _mm512_add_ps(PI,_mm512_mul_ps(cphi1,cphi1));
                      sphi1= xsinf(phi1);
                      sphi1s = _mm512_mul_ps(sphi1,sphi1);
                      trm2 = _mm512_fmadd_ps(k02a2,sphi1s,_mm512_mul_ps(k02b2,cphi1s));
                      x0   = _mm512_pow_ps(trm2,c0);
                      root6= _mm512_rcp14_ps(x0);
                      m    = _mm512_cmp_mask_ps(_mm512_abs_ps(trm1),root6,_CMP_LT_OQ);
                      return (m);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 T_f4423_zmm16r4(const __m512 a,
                                          const __m512 b,
                                          const __m512 phi1,
                                          const __m512 phi2,
                                          const __m512 k0,
                                          bool & status) {

                          using namespace gms::math;
                          __mmask16 m = TM_f4423_helper_zmm16r4(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          register __m512 T,k0c,c,alp,a2,b2;
                          register __m512 sphi,sphi2,cphi,cphi2;
                          register __m512 arg,sarg,rat,x0;
                          a2   = _mm512_mul_ps(a,a);
                          alp  = _mm512_add_ps(PI,_mm512_sub_ps(phi2,phi1));
                          b2   = _mm512_mul_ps(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          x0   = _mm512_fmadd_ps(a2,cphi2,_mm512_mul_ps(b2,sphi2));
                          c    = _mm512_sqrt_ps(x0);
                          k0c  = _mm512_mul_ps(k0,c);
                          arg  = _mm512_mul_ps(k0c,alp);
                          sarg = xsinf(arg);
                          k0c  = negate_zmm16r4(k0c);
                          rat  = _mm512_div_ps(sarg,arg);
                          T    = _mm512_mul_ps(k0c,rat);
                          status = true;
                          return (T);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 T_f4423_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                          const float * __restrict __ATTR_ALIGN__(64) pb,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                          const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                          const float * __restrict __ATTR_ALIGN__(64) pk0,
                                          bool & status) {

                          using namespace gms::math;
                          register __m512 phi1 = _mm512_load_ps(&phi1[0]);
                          register __m512 phi2 = _mm512_load_ps(&phi2[0]);
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm16r4(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          register __m512 T,k0c,c,alp,a2,b2;
                          register __m512 sphi,sphi2,cphi,cphi2;
                          register __m512 arg,sarg,rat,x0;
                          a2   = _mm512_mul_ps(a,a);
                          alp  = _mm512_add_ps(PI,_mm512_sub_ps(phi2,phi1));
                          b2   = _mm512_mul_ps(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          x0   = _mm512_fmadd_ps(a2,cphi2,_mm512_mul_ps(b2,sphi2));
                          c    = _mm512_sqrt_ps(x0);
                          k0c  = _mm512_mul_ps(k0,c);
                          arg  = _mm512_mul_ps(k0c,alp);
                          sarg = xsinf(arg);
                          k0c  = negate_zmm16r4(k0c);
                          rat  = _mm512_div_ps(sarg,arg);
                          T    = _mm512_mul_ps(k0c,rat);
                          status = true;
                          return (T);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 T_f4423_zmm16r4_u(const float * __restrict pa,
                                          const float * __restrict  pb,
                                          const float * __restrict  pphi1,
                                          const float * __restrict  pphi2,
                                          const float * __restrict  pk0,
                                          bool & status) {

                          using namespace gms::math;
                          register __m512 phi1 = _mm512_loadu_ps(&phi1[0]);
                          register __m512 phi2 = _mm512_loadu_ps(&phi2[0]);
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm16r4(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          register __m512 T,k0c,c,alp,a2,b2;
                          register __m512 sphi,sphi2,cphi,cphi2;
                          register __m512 arg,sarg,rat,x0;
                          a2   = _mm512_mul_ps(a,a);
                          alp  = _mm512_add_ps(PI,_mm512_sub_ps(phi2,phi1));
                          b2   = _mm512_mul_ps(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          x0   = _mm512_fmadd_ps(a2,cphi2,_mm512_mul_ps(b2,sphi2));
                          c    = _mm512_sqrt_ps(x0);
                          k0c  = _mm512_mul_ps(k0,c);
                          arg  = _mm512_mul_ps(k0c,alp);
                          sarg = xsinf(arg);
                          k0c  = negate_zmm16r4(k0c);
                          rat  = _mm512_div_ps(sarg,arg);
                          T    = _mm512_mul_ps(k0c,rat);
                          status = true;
                          return (T);
                 }


                   /*
                          Scattering width near the forward direction.
                          Formula 4.4-24

                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4424_zmm16r4(const __m512 a,
                                          const __m512 b,
                                          const __m512 phi1,
                                          const __m512 phi2,
                                          const __m512 k0,
                                          bool & status) {

                         
                          __mmask16 m = TM_f4423_helper_zmm16r4(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,k0c,c,alp,a2,b2;
                          register __m512 sphi,sphi2,cphi,cphi2;
                          register __m512 arg,sarg,rat,x0,x1,x2;
                          a2   = _mm512_mul_ps(a,a);
                          alp  = _mm512_add_ps(PI,_mm512_sub_ps(phi2,phi1));
                          b2   = _mm512_mul_ps(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          x0   = _mm512_fmadd_ps(a2,cphi2,_mm512_mul_ps(b2,sphi2));
                          c    = _mm512_sqrt_ps(x0);
                          k0c  = _mm512_mul_ps(k0,c);
                          arg  = _mm512_mul_ps(k0c,alp);
                          sarg = xsinf(arg);
                          x1   = _mm512_mul_ps(_4,_mm512_mul_ps(k0c,k0c));
                          rat  = _mm512_div_ps(sarg,arg);
                          x2   = _mm512_mul_ps(rat,rat);
                          rcs  = _mm512_mul_ps(k0c,rat);
                          status = true;
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4424_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              bool & status) {

                          
                          register __m512 phi1 = _mm512_load_ps(&phi1[0]);
                          register __m512 phi2 = _mm512_load_ps(&phi2[0]);
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm16r4(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,k0c,c,alp,a2,b2;
                          register __m512 sphi,sphi2,cphi,cphi2;
                          register __m512 arg,sarg,rat,x0,x1,x2;
                          a2   = _mm512_mul_ps(a,a);
                          alp  = _mm512_add_ps(PI,_mm512_sub_ps(phi2,phi1));
                          b2   = _mm512_mul_ps(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          x0   = _mm512_fmadd_ps(a2,cphi2,_mm512_mul_ps(b2,sphi2));
                          c    = _mm512_sqrt_ps(x0);
                          k0c  = _mm512_mul_ps(k0,c);
                          arg  = _mm512_mul_ps(k0c,alp);
                          sarg = xsinf(arg);
                          x1   = _mm512_mul_ps(_4,_mm512_mul_ps(k0c,k0c));
                          rat  = _mm512_div_ps(sarg,arg);
                          x2   = _mm512_mul_ps(rat,rat);
                          rcs  = _mm512_mul_ps(k0c,rat);
                          status = true;
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4424_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pphi1,
                                              const float * __restrict  pphi2,
                                              const float * __restrict  pk0,
                                              bool & status) {

                          
                          register __m512 phi1 = _mm512_loadu_ps(&phi1[0]);
                          register __m512 phi2 = _mm512_loadu_ps(&phi2[0]);
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm16r4(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,k0c,c,alp,a2,b2;
                          register __m512 sphi,sphi2,cphi,cphi2;
                          register __m512 arg,sarg,rat,x0,x1,x2;
                          a2   = _mm512_mul_ps(a,a);
                          alp  = _mm512_add_ps(PI,_mm512_sub_ps(phi2,phi1));
                          b2   = _mm512_mul_ps(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          x0   = _mm512_fmadd_ps(a2,cphi2,_mm512_mul_ps(b2,sphi2));
                          c    = _mm512_sqrt_ps(x0);
                          k0c  = _mm512_mul_ps(k0,c);
                          arg  = _mm512_mul_ps(k0c,alp);
                          sarg = xsinf(arg);
                          x1   = _mm512_mul_ps(_4,_mm512_mul_ps(k0c,k0c));
                          rat  = _mm512_div_ps(sarg,arg);
                          x2   = _mm512_mul_ps(rat,rat);
                          rcs  = _mm512_mul_ps(k0c,rat);
                          status = true;
                          return (rcs);
                 }


                   /*
                         Scattering width in the exact forward direction (alpha == 0).
                         Formula 4.4-25
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4425_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                            const __m512 phi) {

                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,_4k0,a2,b2,sphi,sphi2;
                          register __m512 cphi,cphi2;
                          a2   = _mm512_mul_ps(a,a);
                          sphi = xsinf(phi);
                          b2   = _mm512_mul_ps(b,b);
                          cphi = xcosf(phi);
                          _4k0 = _mm512_mul_ps(_4,k0);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          x0   = _mm512_fmadd_ps(a2,sphi2,_mm512_mul_ps(b2,cphi2));
                          rcs  = _mm512_mul_ps(_4k0,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4425_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const float * __restrict __ATTR_ALIGN__(64)  pa,
                                              const float * __restrict __ATTR_ALIGN__(64)  pb,
                                              const float * __restrict __ATTR_ALIGN__(64)  pphi) {

                          register __m512 phi2 = _mm512_load_ps(&phi2[0]);
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 b    = _mm512_load_ps(&pb[0]);
                          register __m512 k0   = _mm512_load_ps(&pk0[0]); 
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,_4k0,a2,b2,sphi,sphi2;
                          register __m512 cphi,cphi2;
                          a2   = _mm512_mul_ps(a,a);
                          sphi = xsinf(phi);
                          b2   = _mm512_mul_ps(b,b);
                          cphi = xcosf(phi);
                          _4k0 = _mm512_mul_ps(_4,k0);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          x0   = _mm512_fmadd_ps(a2,sphi2,_mm512_mul_ps(b2,cphi2));
                          rcs  = _mm512_mul_ps(_4k0,x0);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4425_zmm16r4_u(const float * __restrict   pk0,
                                              const float * __restrict   pa,
                                              const float * __restrict   pb,
                                              const float * __restrict  pphi) {

                          register __m512 phi2 = _mm512_loadu_ps(&phi2[0]);
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 b    = _mm512_loadu_ps(&pb[0]);
                          register __m512 k0   = _mm512_loadu_ps(&pk0[0]); 
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs,_4k0,a2,b2,sphi,sphi2;
                          register __m512 cphi,cphi2;
                          a2   = _mm512_mul_ps(a,a);
                          sphi = xsinf(phi);
                          b2   = _mm512_mul_ps(b,b);
                          cphi = xcosf(phi);
                          _4k0 = _mm512_mul_ps(_4,k0);
                          cphi2= _mm512_mul_ps(cphi,cphi);
                          sphi2= _mm512_mul_ps(sphi,sphi);
                          x0   = _mm512_fmadd_ps(a2,sphi2,_mm512_mul_ps(b2,cphi2));
                          rcs  = _mm512_mul_ps(_4k0,x0);
                          return (rcs);
                 }


                    /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          TM-case, formula 4.4-26
                     */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4426_zmm16r4(const __m512 k0,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 phi1,
                                         const __m512 phi2,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui,
                                         __m512 * __restrict TMr,
                                         __m512 * __restrict TMi) {

                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,murpba,muipba;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512 facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        faci   = _mm512_mul_ps(pi4,_mm512_mul_ps(k0a2,ba));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        *TMr   = _mm512_mul_ps(facr,tmpr);
                        *TMi   = _mm512_mul_ps(faci,tmpi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4426_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pb,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                         const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                         const float * __restrict __ATTR_ALIGN__(64) pmur,
                                         const float * __restrict __ATTR_ALIGN__(64) pmui,
                                         float * __restrict __ATTR_ALIGN__(64) TMr,
                                         float * __restrict __ATTR_ALIGN__(64) TMi) {

                        register __m512 k0    = _mm512_load_ps(&pk0[0]);
                        register __m512 a     = _mm512_load_ps(&pa[0]);
                        register __m512 b     = _mm512_load_ps(&pb[0]);
                        register __m512 phi1  = _mm512_load_ps(&pphi1[0]);
                        register __m512 phi2  = _mm512_load_ps(&pphi2[0]);
                        register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_load_ps(&pmur[0]);
                        register __m512 mui   = _mm512_load_ps(&pmui[0]);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,murpba,muipba;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512 facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        faci   = _mm512_mul_ps(pi4,_mm512_mul_ps(k0a2,ba));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_store_ps(&TMr[0]  ,_mm512_mul_ps(facr,tmpr));
                        _mm512_store_ps(&TMi[0]  ,_mm512_mul_ps(faci,tmpi));
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4426_zmm16r4_u(const float * __restrict  pk0,
                                         const float * __restrict  pa,
                                         const float * __restrict  pb,
                                         const float * __restrict  pphi1,
                                         const float * __restrict pphi2,
                                         const float * __restrict  pepsr,
                                         const float * __restrict  pepsi,
                                         const float * __restrict  pmur,
                                         const float * __restrict  pmui,
                                         float * __restrict  TMr,
                                         float * __restrict  TMi) {

                        register __m512 k0    = _mm512_loadu_ps(&pk0[0]);
                        register __m512 a     = _mm512_loadu_ps(&pa[0]);
                        register __m512 b     = _mm512_loadu_ps(&pb[0]);
                        register __m512 phi1  = _mm512_loadu_ps(&pphi1[0]);
                        register __m512 phi2  = _mm512_loadu_ps(&pphi2[0]);
                        register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui   = _mm512_loadu_ps(&pmui[0]);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,murpba,muipba;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512 facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        faci   = _mm512_mul_ps(pi4,_mm512_mul_ps(k0a2,ba));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_storeu_ps(&TMr[0]  ,_mm512_mul_ps(facr,tmpr));
                        _mm512_storeu_ps(&TMi[0]  ,_mm512_mul_ps(faci,tmpi));
                }


                   /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          TE-case, formula 4.4-27
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4427_zmm16r4(const __m512 k0,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 phi1,
                                         const __m512 phi2,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui,
                                         __m512 * __restrict TEr,
                                         __m512 * __restrict TEi) {

                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512 facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        faci   = _mm512_mul_ps(pi4,_mm512_mul_ps(k0a2,ba));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        *TEr   = _mm512_mul_ps(facr,tmpr);
                        *TEi   = _mm512_mul_ps(faci,tmpi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4427_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pb,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                         const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                         const float * __restrict __ATTR_ALIGN__(64) pmur,
                                         const float * __restrict __ATTR_ALIGN__(64) pmui,
                                         float * __restrict __ATTR_ALIGN__(64) TEr,
                                         float * __restrict __ATTR_ALIGN__(64) TEi) {

                        register __m512 k0    = _mm512_load_ps(&pk0[0]);
                        register __m512 a     = _mm512_load_ps(&pa[0]);
                        register __m512 b     = _mm512_load_ps(&pb[0]);
                        register __m512 phi1  = _mm512_load_ps(&pphi1[0]);
                        register __m512 phi2  = _mm512_load_ps(&pphi2[0]);
                        register __m512 epsr  = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_load_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_load_ps(&pmur[0]);
                        register __m512 mui   = _mm512_load_ps(&pmui[0]);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512 facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        faci   = _mm512_mul_ps(pi4,_mm512_mul_ps(k0a2,ba));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_store_ps(&TEr[0] ,_mm512_mul_ps(facr,tmpr));
                        _mm512_store_ps(&TEi[0] ,_mm512_mul_ps(faci,tmpi));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4427_zmm16r4_u(const float * __restrict  pk0,
                                         const float * __restrict  pa,
                                         const float * __restrict  pb,
                                         const float * __restrict  pphi1,
                                         const float * __restrict pphi2,
                                         const float * __restrict  pepsr,
                                         const float * __restrict  pepsi,
                                         const float * __restrict pmur,
                                         const float * __restrict  pmui,
                                         float * __restrict  TEr,
                                         float * __restrict  TEi) {

                        register __m512 k0    = _mm512_loadu_ps(&pk0[0]);
                        register __m512 a     = _mm512_loadu_ps(&pa[0]);
                        register __m512 b     = _mm512_loadu_ps(&pb[0]);
                        register __m512 phi1  = _mm512_loadu_ps(&pphi1[0]);
                        register __m512 phi2  = _mm512_loadu_ps(&pphi2[0]);
                        register __m512 epsr  = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi  = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 mur   = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui   = _mm512_loadu_ps(&pmui[0]);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512 facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        faci   = _mm512_mul_ps(pi4,_mm512_mul_ps(k0a2,ba));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_storeu_ps(&TEr[0] ,_mm512_mul_ps(facr,tmpr));
                        _mm512_storeu_ps(&TEi[0] ,_mm512_mul_ps(faci,tmpi));
                }


                  /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Bistatic scattering width (RCS).
                          TM-case.
                          Formula 4.4-28
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4428_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                            const __m512 phi1,
                                            const __m512 phi2,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {
                                        
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512 cphi2,cphi1,sphi2,sphi1,murpba,muipba,b2,a2,pia;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512 fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_ps(b,b);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        pia    = _mm512_mul_ps(PI,a);
                        a2     = _mm512_mul_ps(a,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_ps(b2,a2)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        k0a3   = _mm512_mul_ps(k0a2,k0a);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4428_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {
                                        
                        register __m512  k0   = _mm512_load_ps(&pk0[0]);
                        register __m512  a    = _mm512_load_ps(&pa[0]);
                        register __m512  b    = _mm512_load_ps(&pb[0]);
                        register __m512  phi1 = _mm512_load_ps(&pphi1[0]);
                        register __m512  phi2 = _mm512_load_ps(&pphi2[0]);
                        register __m512  epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_load_ps(&pmur[0]);
                        register __m512  mui  = _mm512_load_ps(&pmui[0]);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512 cphi2,cphi1,sphi2,sphi1,murpba,muipba,b2,a2,pia;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512 fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_ps(b,b);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        pia    = _mm512_mul_ps(PI,a);
                        a2     = _mm512_mul_ps(a,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_ps(b2,a2)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        k0a3   = _mm512_mul_ps(k0a2,k0a);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4428_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pphi1,
                                              const float * __restrict  pphi2,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {
                                        
                        register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512  a    = _mm512_loadu_ps(&pa[0]);
                        register __m512  b    = _mm512_loadu_ps(&pb[0]);
                        register __m512  phi1 = _mm512_loadu_ps(&pphi1[0]);
                        register __m512  phi2 = _mm512_loadu_ps(&pphi2[0]);
                        register __m512  epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512  mui  = _mm512_loadu_ps(&pmui[0]);
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512 cphi2,cphi1,sphi2,sphi1,murpba,muipba,b2,a2,pia;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512 fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_ps(b,b);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        pia    = _mm512_mul_ps(PI,a);
                        a2     = _mm512_mul_ps(a,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_ps(b2,a2)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        k0a3   = _mm512_mul_ps(k0a2,k0a);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                  /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Bistatic scattering width (RCS).
                          TE-case.
                          Formula 4.4-29

                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4429_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                         const __m512 phi1,
                                         const __m512 phi2,
                                         const __m512 epsr,
                                         const __m512 epsi,
                                         const __m512 mur,
                                         const __m512 mui) {
                                         
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba,cabs;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512 fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_ps(a,a);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_ps(b,b);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        b2a2   = _mm512_div_ps(b2,a2);
                        pia    = _mm512_mul_ps(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        k0a3   = _mm512_mul_ps(k0a,k0a2);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4429_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {
                              
                        register __m512  k0   = _mm512_load_ps(&pk0[0]);
                        register __m512  a    = _mm512_load_ps(&pa[0]);
                        register __m512  b    = _mm512_load_ps(&pb[0]);
                        register __m512  phi1 = _mm512_load_ps(&pphi1[0]);
                        register __m512  phi2 = _mm512_load_ps(&pphi2[0]);
                        register __m512  epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_load_ps(&pmur[0]);
                        register __m512  mui  = _mm512_load_ps(&pmui[0]);           
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba,cabs;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512 fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_ps(a,a);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_ps(b,b);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        b2a2   = _mm512_div_ps(b2,a2);
                        pia    = _mm512_mul_ps(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        k0a3   = _mm512_mul_ps(k0a,k0a2);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4429_zmm16r4_u(const float * __restrict pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict pb,
                                              const float * __restrict pphi1,
                                              const float * __restrict  pphi2,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {
                              
                        register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512  a    = _mm512_loadu_ps(&pa[0]);
                        register __m512  b    = _mm512_loadu_ps(&pb[0]);
                        register __m512  phi1 = _mm512_loadu_ps(&pphi1[0]);
                        register __m512  phi2 = _mm512_loadu_ps(&pphi2[0]);
                        register __m512  epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512  mui  = _mm512_loadu_ps(&pmui[0]);           
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba,cabs;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512 fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_ps(a,a);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_ps(b,b);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        b2a2   = _mm512_div_ps(b2,a2);
                        pia    = _mm512_mul_ps(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        k0a3   = _mm512_mul_ps(k0a,k0a2);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_ps(cphi2,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        sphit  = _mm512_mul_ps(sphi2,sphi1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphit,murpba);
                        t1i    = _mm512_div_ps(cphit,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphit,murmba);
                        t2i    = _mm512_div_ps(sphit,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                    /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Backscattering  width (RCS).
                          TM-case.
                          Formula 4.4-30
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4430_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                            const __m512 phi1,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {
                                        
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512 cphi1,sphi1,cphi1s,sphi1s,murpba,muipba,b2,a2,pia;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512 fac,_1ba,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_ps(b,b);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        pia    = _mm512_mul_ps(PI,a);
                        a2     = _mm512_mul_ps(a,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_ps(b2,a2)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi1s = _mm512_mul_ps(cphi1,cphi1);
                        k0a3   = _mm512_mul_ps(k0a2,k0a);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi1s  = _mm512_mul_ps(sphi1,sphi1);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphi1s,murpba);
                        t1i    = _mm512_div_ps(cphi1s,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphi1s,murmba);
                        t2i    = _mm512_div_ps(sphi1s,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4430_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {
                                 
                                 
                        register __m512  k0   = _mm512_load_ps(&pk0[0]);
                        register __m512  a    = _mm512_load_ps(&pa[0]);
                        register __m512  b    = _mm512_load_ps(&pb[0]);
                        register __m512  phi1 = _mm512_load_ps(&pphi1[0]);
                        register __m512  epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_load_ps(&pmur[0]);
                        register __m512  mui  = _mm512_load_ps(&pmui[0]);              
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512 cphi1,sphi1,cphi1s,sphi1s,murpba,muipba,b2,a2,pia;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512 fac,_1ba,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_ps(b,b);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        pia    = _mm512_mul_ps(PI,a);
                        a2     = _mm512_mul_ps(a,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_ps(b2,a2)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi1s = _mm512_mul_ps(cphi1,cphi1);
                        k0a3   = _mm512_mul_ps(k0a2,k0a);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi1s  = _mm512_mul_ps(sphi1,sphi1);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphi1s,murpba);
                        t1i    = _mm512_div_ps(cphi1s,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphi1s,murmba);
                        t2i    = _mm512_div_ps(sphi1s,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4430_zmm16r4_u(const float * __restrict pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict pphi1,
                                              const float * __restrict pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {
                                 
                                 
                        register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512  a    = _mm512_loadu_ps(&pa[0]);
                        register __m512  b    = _mm512_loadu_ps(&pb[0]);
                        register __m512  phi1 = _mm512_loadu_ps(&pphi1[0]);
                        register __m512  epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512  mui  = _mm512_loadu_ps(&pmui[0]);              
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512 cphi1,sphi1,cphi1s,sphi1s,murpba,muipba,b2,a2,pia;
                        register __m512 murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512 fac,_1ba,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_ps(b,b);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_ps(b,a);
                        pia    = _mm512_mul_ps(PI,a);
                        a2     = _mm512_mul_ps(a,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_ps(b2,a2)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        cphi1s = _mm512_mul_ps(cphi1,cphi1);
                        k0a3   = _mm512_mul_ps(k0a2,k0a);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(epsrm1,murm1);
                        sphi1s  = _mm512_mul_ps(sphi1,sphi1);
                        t0i    = _mm512_sub_ps(epsim1,muim1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        murpba = _mm512_add_ps(mur,ba);
                        muipba = _mm512_add_ps(mui,ba);
                        t1r    = _mm512_div_ps(cphi1s,murpba);
                        t1i    = _mm512_div_ps(cphi1s,muipba);
                        murmba = _mm512_fmadd_ps(mur,ba,_1);
                        muimba = _mm512_fmadd_ps(mui,ba,_1);
                        t2r    = _mm512_div_ps(sphi1s,murmba);
                        t2i    = _mm512_div_ps(sphi1s,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                  /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Backscattering  width (RCS).
                          TE-case.
                          Formula 4.4-31
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4431_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                            const __m512 phi1,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {
                                         
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi1,sphi1,cphi1s,sphi1s,epsrpba,epsipba,cabs;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512 fac,_1ba,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_ps(a,a);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_ps(b,b);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        b2a2   = _mm512_div_ps(b2,a2);
                        pia    = _mm512_mul_ps(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        k0a3   = _mm512_mul_ps(k0a,k0a2);
                        cphi1s = _mm512_mul_ps(cphi1,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi1s = _mm512_mul_ps(sphi1,sphi1);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphi1s,murpba);
                        t1i    = _mm512_div_ps(cphi1s,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphi1s,murmba);
                        t2i    = _mm512_div_ps(sphi1s,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4431_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {
                                
                        register __m512  k0   = _mm512_load_ps(&pk0[0]);
                        register __m512  a    = _mm512_load_ps(&pa[0]);
                        register __m512  b    = _mm512_load_ps(&pb[0]);
                        register __m512  phi1 = _mm512_load_ps(&pphi1[0]);
                        register __m512  epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_load_ps(&pmur[0]);
                        register __m512  mui  = _mm512_load_ps(&pmui[0]);           
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi1,sphi1,cphi1s,sphi1s,epsrpba,epsipba,cabs;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512 fac,_1ba,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_ps(a,a);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_ps(b,b);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        b2a2   = _mm512_div_ps(b2,a2);
                        pia    = _mm512_mul_ps(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        k0a3   = _mm512_mul_ps(k0a,k0a2);
                        cphi1s = _mm512_mul_ps(cphi1,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi1s = _mm512_mul_ps(sphi1,sphi1);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphi1s,murpba);
                        t1i    = _mm512_div_ps(cphi1s,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphi1s,murmba);
                        t2i    = _mm512_div_ps(sphi1s,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4431_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pphi1,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {
                                
                        register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512  a    = _mm512_loadu_ps(&pa[0]);
                        register __m512  b    = _mm512_loadu_ps(&pb[0]);
                        register __m512  phi1 = _mm512_loadu_ps(&pphi1[0]);
                        register __m512  epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512  epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512  mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512  mui  = _mm512_loadu_ps(&pmui[0]);           
                        const __m512 _1  = _mm512_set1_ps(1.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        register __m512 rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512 cphi1,sphi1,cphi1s,sphi1s,epsrpba,epsipba,cabs;
                        register __m512 epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512 fac,_1ba,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_ps(a,a);
                        k0a    = _mm512_mul_ps(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_ps(b,b);
                        ba     = _mm512_div_ps(b,a);
                        epsrm1 = _mm512_sub_ps(epsr,_1);
                        b2a2   = _mm512_div_ps(b2,a2);
                        pia    = _mm512_mul_ps(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_ps(k0a,k0a);
                        epsim1 = _mm512_sub_ps(epsi,_1);
                        k0a3   = _mm512_mul_ps(k0a,k0a2);
                        cphi1s = _mm512_mul_ps(cphi1,cphi1);
                        murm1  = _mm512_sub_ps(mur,_1);
                        muim1  = _mm512_sub_ps(mui,_1);
                        _1ba   = _mm512_add_ps(_1,ba);
                        t0r    = _mm512_sub_ps(murm1,epsrm1);
                        sphi1s = _mm512_mul_ps(sphi1,sphi1);
                        t0i    = _mm512_sub_ps(muim1,epsim1);
                        fac    = _mm512_mul_ps(_mm512_mul_ps(pia,pi4),
                                               _mm512_mul_ps(k0a3,b2a2));
                        epsrpba= _mm512_add_ps(epsr,ba);
                        epsipba= _mm512_add_ps(epsi,ba);
                        t1r    = _mm512_div_ps(cphi1s,murpba);
                        t1i    = _mm512_div_ps(cphi1s,muipba);
                        epsrmba= _mm512_fmadd_ps(epsr,ba,_1);
                        epsimba= _mm512_fmadd_ps(epsi,ba,_1);
                        t2r    = _mm512_div_ps(sphi1s,murmba);
                        t2i    = _mm512_div_ps(sphi1s,muimba);
                        t3r    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1r,t2r));
                        t3i    = _mm512_mul_ps(_1ba,_mm512_add_ps(t1i,t2i));
                        cmul_zmm16r4(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm16r4(tmpr,tmpi);
                        rcs    = _mm512_mul_ps(fac,cabs);
                        return (rcs);
                }


                 /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Forward scattering (phi2 = pi+phi1)  width (RCS).
                          TM-case.
                          Formula 4.4-32
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4432_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                            const __m512 phi1,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {

                           return (rcs_f4430_zmm16r4(k0,a,b,phi1,epsr,epsi,mur,mui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4432_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {

                           return (rcs_f4430_zmm16r4_a(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4432_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pphi1,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {

                           return (rcs_f4430_zmm16r4_u(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               } 


                 /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Forward scattering (phi2 = pi+phi1)  width (RCS).
                          TE-case.
                          Formula 4.4-33

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4433_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 b,
                                            const __m512 phi1,
                                            const __m512 epsr,
                                            const __m512 epsi,
                                            const __m512 mur,
                                            const __m512 mui) {

                           return (rcs_f4431_zmm16r4(k0,a,b,phi1,epsr,epsi,mur,mui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4433_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const float * __restrict __ATTR_ALIGN__(64) pmur,
                                              const float * __restrict __ATTR_ALIGN__(64) pmui) {

                           return (rcs_f4431_zmm16r4_a(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4433_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pphi1,
                                              const float * __restrict  pepsr,
                                              const float * __restrict  pepsi,
                                              const float * __restrict  pmur,
                                              const float * __restrict  pmui) {

                           return (rcs_f4431_zmm16r4_u(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               }


                








                  

 



                   






                   

 




      } // radiolocation


} // gms









#endif /*__GMS_RCS_CYLINDER_ZMM16R4_HPP__*/
