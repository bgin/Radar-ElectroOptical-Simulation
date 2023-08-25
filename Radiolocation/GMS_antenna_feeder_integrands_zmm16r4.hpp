

#ifndef __GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_HPP__
#define __GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_HPP__


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

    const unsigned int GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_FULLVER =
      1000U*GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_MAJOR+
      100U*GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_MINOR+
      10U*GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_MICRO;
    const char * const GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_CREATION_DATE = "23-08-2023 15:31 PM +00200 (WED 23 AUG 2023 GMT+2)";
    const char * const GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4_DESCRIPTION   = "AVX512 (single) optimized antenna and feeder integrands.";

}


#include <immintrin.h>
#include <complex>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"
#include "GMS_em_fields_zmm16r4.hpp"
#include "GMS_cephes.h"

namespace gms {

      
             namespace radiolocation {
       
#if !defined(__ANTENNA_FEEDER_PF_CACHE_HINT__)
#define  __ANTENNA_FEEDER_PF_CACHE_HINT__ 1
#endif                  
                 
                 struct __ATTR_ALIGN__(64) fwork_t {
                   
                          float * __restrict pxr;
                          float * __restrict pxi;
                          float * __restrict pyr;
                          float * __restrict pyi;
                          float * __restrict pzr;
                          float * __restrict pzi;
                          PAD_TO(0,16)
                   };
                   
                   
                  // Spherical unit vector.
                  struct __ATTR_ALIGN__(64) SUV_zmm16r4_t {
                  
                         __m512 x;
                         __m512 y;
                         __m512 z;
                  };
                  
                  
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           void  f2135_integrand_zmm16r4_u6x_a(const float * __restrict __ATTR_ALIGN__(64) pxre,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pxim,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pyre,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pyim,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pzre,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pzim,
	                                              fwork_t fw, //work arrays (caller allocated)
	                                              const __m512 cer,
	                                              const __m512 cei,
	                                              const std::complex<float> tmp,
                                                      const int32_t n,
                                                      const int32_t  PF_DIST) {
                                                      
                           if(__builtin_expect(n<=0,0)) { return;}
                           register __m512 xre,xim,yre,yim,zre,zim;
                            __m512 t0r,t0i,t1r,t1i,t2r,t2i;
                           int32_t i;
                           
                           for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_NTA);
#endif
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_load_ps(&pxre[i+16]);
                             xim = _mm512_load_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_load_ps(&pyre[i+16]);
                             yim = _mm512_load_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_load_ps(&pzre[i+16]);
                             zim = _mm512_load_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+16],t2r);
                             _mm512_store_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_load_ps(&pxre[i+32]);
                             xim = _mm512_load_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+32],t0r);
                             _mm512_store_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_load_ps(&pyre[i+32]);
                             yim = _mm512_load_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+32],t1r);
                             _mm512_store_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_load_ps(&pzre[i+32]);
                             zim = _mm512_load_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+32],t2r);
                             _mm512_store_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_load_ps(&pxre[i+48]);
                             xim = _mm512_load_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+48],t0r);
                             _mm512_store_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_load_ps(&pyre[i+48]);
                             yim = _mm512_load_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+48],t1r);
                             _mm512_store_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_load_ps(&pzre[i+48]);
                             zim = _mm512_load_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+48],t2r);
                             _mm512_store_ps(&fw.zyi[i+48],t2i)
                             xre = _mm512_load_ps(&pxre[i+64]);
                             xim = _mm512_load_ps(&pxim[i+64]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+64],t0r);
                             _mm512_store_ps(&fw.pxi[i+64],t0i);
                             yre = _mm512_load_ps(&pyre[i+64]);
                             yim = _mm512_load_ps(&pyim[i+64]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+64],t1r);
                             _mm512_store_ps(&fw.pyi[i+64],t1i);
                             zre = _mm512_load_ps(&pzre[i+64]);
                             zim = _mm512_load_ps(&pzim[i+64]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+64],t2r);
                             _mm512_store_ps(&fw.zyi[i+64],t2i)
                             xre = _mm512_load_ps(&pxre[i+80]);
                             xim = _mm512_load_ps(&pxim[i+80]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+80],t0r);
                             _mm512_store_ps(&fw.pxi[i+80],t0i);
                             yre = _mm512_load_ps(&pyre[i+80]);
                             yim = _mm512_load_ps(&pyim[i+80]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+80],t1r);
                             _mm512_store_ps(&fw.pyi[i+80],t1i);
                             zre = _mm512_load_ps(&pzre[i+80]);
                             zim = _mm512_load_ps(&pzim[i+80]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+80],t2r);
                             _mm512_store_ps(&fw.zyi[i+80],t2i)
                        }  
                        
                        for(; (i+63) < n; i += 64) {
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_load_ps(&pxre[i+16]);
                             xim = _mm512_load_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_load_ps(&pyre[i+16]);
                             yim = _mm512_load_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_load_ps(&pzre[i+16]);
                             zim = _mm512_load_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+16],t2r);
                             _mm512_store_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_load_ps(&pxre[i+32]);
                             xim = _mm512_load_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+32],t0r);
                             _mm512_store_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_load_ps(&pyre[i+32]);
                             yim = _mm512_load_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+32],t1r);
                             _mm512_store_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_load_ps(&pzre[i+32]);
                             zim = _mm512_load_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+32],t2r);
                             _mm512_store_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_load_ps(&pxre[i+48]);
                             xim = _mm512_load_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+48],t0r);
                             _mm512_store_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_load_ps(&pyre[i+48]);
                             yim = _mm512_load_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+48],t1r);
                             _mm512_store_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_load_ps(&pzre[i+48]);
                             zim = _mm512_load_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+48],t2r);
                             _mm512_store_ps(&fw.zyi[i+48],t2i)
                      }   
                      
                      for(; (i+31) < n; i += 32) {
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_load_ps(&pxre[i+16]);
                             xim = _mm512_load_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_load_ps(&pyre[i+16]);
                             yim = _mm512_load_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_load_ps(&pzre[i+16]);
                             zim = _mm512_load_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+16],t2r);
                             _mm512_store_ps(&fw.zyi[i+16],t2i) 
                      }   
                      
                      for(; (i+15) < n; i += 16) {
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                      }  
                      
                     
                      for(; (i+0) < n; i += 1) {
                           register float xr               = pxre[i];
                           register float xi               = pxim[i];
                           register std::complex<float> cx = {xr,xi};
                           register std::complex<float> xx = cx*tmp2;
                           fw.pxr[i]                       = xx.real(); 
                           fw.pxi[i]                       = xx.imag();
                           register float yr               = pyre[i];
                           register float yi               = pyim[i];
                           register std::complex<float> cy = {yr,yi};
                           register std::complex<float> yy = cy*tmp2;
                           fw.pyr[i]                       = yy.real(); 
                           fw.pyi[i]                       = yy.imag();
                           register float zr               = pzre[i];
                           register float zi               = pzim[i];
                           register std::complex<float> cz = {zr,zi};
                           register std::complex<float> zz = cz*tmp2;
                           fw.pzr[i]                       = zz.real(); 
                           fw.pzi[i]                       = zz.imag();
                      }  
                      
                                                       
                  }
                  
                  
                      
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           void  f2135_integrand_zmm16r4_u6x_u(const float * __restrict  pxre,
	                                              const float * __restrict  pxim,
	                                              const float * __restrict  pyre,
	                                              const float * __restrict  pyim,
	                                              const float * __restrict  pzre,
	                                              const float * __restrict  pzim,
	                                              fwork_t fw, //work arrays (caller allocated)
	                                              const __m512 cer,
	                                              const __m512 cei,
	                                              const std::complex<float> tmp,
                                                      const int32_t n,
                                                      const int32_t  PF_DIST) {
                                                      
                                                      
                           if(__builtin_expect(n<=0,0)) { return;}
                           register __m512 xre,xim,yre,yim,zre,zim;
                            __m512 t0r,t0i,t1r,t1i,t2r,t2i;
                           int32_t i;
                           
                           for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_NTA);
#endif
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+16]);
                             xim = _mm512_loadu_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+16]);
                             yim = _mm512_loadu_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+16]);
                             zim = _mm512_loadu_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+16],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+32]);
                             xim = _mm512_loadu_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+32]);
                             yim = _mm512_loadu_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+32]);
                             zim = _mm512_loadu_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+32],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+48]);
                             xim = _mm512_loadu_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+48]);
                             yim = _mm512_loadu_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+48]);
                             zim = _mm512_loadu_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+48],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+48],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+64]);
                             xim = _mm512_loadu_ps(&pxim[i+64]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+64],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+64],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+64]);
                             yim = _mm512_loadu_ps(&pyim[i+64]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+64],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+64],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+64]);
                             zim = _mm512_loadu_ps(&pzim[i+64]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+64],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+64],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+80]);
                             xim = _mm512_loadu_ps(&pxim[i+80]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+80],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+80],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+80]);
                             yim = _mm512_loadu_ps(&pyim[i+80]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+80],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+80],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+80]);
                             zim = _mm512_loadu_ps(&pzim[i+80]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+80],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+80],t2i)
                        }  
                        
                        for(; (i+63) < n; i += 64) {
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+16]);
                             xim = _mm512_loadu_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+16]);
                             yim = _mm512_loadu_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+16]);
                             zim = _mm512_loadu_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+16],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+32]);
                             xim = _mm512_loadu_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+32]);
                             yim = _mm512_loadu_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+32]);
                             zim = _mm512_loadu_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+32],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+48]);
                             xim = _mm512_loadu_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+48]);
                             yim = _mm512_loadu_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+48]);
                             zim = _mm512_loadu_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+48],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+48],t2i)
                      }   
                      
                      for(; (i+31) < n; i += 32) {
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+16]);
                             xim = _mm512_loadu_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+16]);
                             yim = _mm512_loadu_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+16]);
                             zim = _mm512_loadu_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+16],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+16],t2i) 
                      }   
                      
                      for(; (i+15) < n; i += 16) {
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);  
                      }  
                      
                     
                      for(; (i+0) < n; i += 1) {
                           register float xr               = pxre[i];
                           register float xi               = pxim[i];
                           register std::complex<float> cx = {xr,xi};
                           register std::complex<float> xx = cx*tmp2;
                           fw.pxr[i]                       = xx.real(); 
                           fw.pxi[i]                       = xx.imag();
                           register float yr               = pyre[i];
                           register float yi               = pyim[i];
                           register std::complex<float> cy = {yr,yi};
                           register std::complex<float> yy = cy*tmp2;
                           fw.pyr[i]                       = yy.real(); 
                           fw.pyi[i]                       = yy.imag();
                           register float zr               = pzre[i];
                           register float zi               = pzim[i];
                           register std::complex<float> cz = {zr,zi};
                           register std::complex<float> zz = cz*tmp2;
                           fw.pzr[i]                       = zz.real(); 
                           fw.pzi[i]                       = zz.imag();
                      }  
                      
                                                       
                  }
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           void f2235_integrand_zmm16r4_u6x_a(const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                              const float * __restrict __ATTR_ALIGN__(64) prho,
	                                              const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                              fwork_t fw, ////work arrays (caller allocated)
	                                              const float k,
	                                              const int32_t n,
	                                              const int32_t PF_DIST) {
	                                              
	                 if(__builtin_expect(n<=0,0)) { return;}
	                 register __m512 jxr,jxi,jyr,jyi,jzr,jzi;
	                 register __m512 rho,cst,ir,ii;
	                 register __m512 ear,eai,cer,cei;
	                 register __m512 vk,t0;
	                 __m512   t0r,t0i,t1r,t1i,t2r,t2i;
	                 int32_t i;
	                 
	                 vk = _mm512_set1_ps(k);
	                 ii = _mm512_set1_ps(1.0f);
	                 ir = _mm512_setzero_ps();
	                 ear= ir;
	                 t0 = _mm512_mul_ps(ii,vk); // shall be removed.
	                 for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_NTA);
#endif	                 
                              cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_load_ps(&pcst[i+16]);
                              rho = _mm512_load_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+16]);
                              jxi = _mm512_load_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+16],t0r);
                              _mm512_store_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+16]);
                              jyi = _mm512_load_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+16],t1r);
                              _mm512_store_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+16]);
                              jzi = _mm512_load_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+16],t2r);
                              _mm512_store_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_load_ps(&pcst[i+32]);
                              rho = _mm512_load_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+32]);
                              jxi = _mm512_load_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+32],t0r);
                              _mm512_store_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+32]);
                              jyi = _mm512_load_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+32],t1r);
                              _mm512_store_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+32]);
                              jzi = _mm512_load_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+32],t2r);
                              _mm512_store_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_load_ps(&pcst[i+48]);
                              rho = _mm512_load_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+48]);
                              jxi = _mm512_load_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+48],t0r);
                              _mm512_store_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+48]);
                              jyi = _mm512_load_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+48],t1r);
                              _mm512_store_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+48]);
                              jzi = _mm512_load_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+48],t2r);
                              _mm512_store_ps(&fw.pzi[i+48],t2i);
                              cst = _mm512_load_ps(&pcst[i+64]);
                              rho = _mm512_load_ps(&prho[i+64]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+64]);
                              jxi = _mm512_load_ps(&pjxi[i+64]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+64],t0r);
                              _mm512_store_ps(&fw.pxi[i+64],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+64]);
                              jyi = _mm512_load_ps(&pjyi[i+64]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+64],t1r);
                              _mm512_store_ps(&fw.pyi[i+64],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+64]);
                              jzi = _mm512_load_ps(&pjzi[i+64]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+64],t2r);
                              _mm512_store_ps(&fw.pzi[i+64],t2i);
                              cst = _mm512_load_ps(&pcst[i+80]);
                              rho = _mm512_load_ps(&prho[i+80]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+80]);
                              jxi = _mm512_load_ps(&pjxi[i+80]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+80],t0r);
                              _mm512_store_ps(&fw.pxi[i+80],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+80]);
                              jyi = _mm512_load_ps(&pjyi[i+80]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+80],t1r);
                              _mm512_store_ps(&fw.pyi[i+80],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+80]);
                              jzi = _mm512_load_ps(&pjzi[i+80]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+80],t2r);
                              _mm512_store_ps(&fw.pzi[i+80],t2i);
	                 }  
	                 
	                 for(; (i+63) < n; i += 64) {
	                      cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_load_ps(&pcst[i+16]);
                              rho = _mm512_load_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+16]);
                              jxi = _mm512_load_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+16],t0r);
                              _mm512_store_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+16]);
                              jyi = _mm512_load_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+16],t1r);
                              _mm512_store_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+16]);
                              jzi = _mm512_load_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+16],t2r);
                              _mm512_store_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_load_ps(&pcst[i+32]);
                              rho = _mm512_load_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+32]);
                              jxi = _mm512_load_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+32],t0r);
                              _mm512_store_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+32]);
                              jyi = _mm512_load_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+32],t1r);
                              _mm512_store_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+32]);
                              jzi = _mm512_load_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+32],t2r);
                              _mm512_store_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_load_ps(&pcst[i+48]);
                              rho = _mm512_load_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+48]);
                              jxi = _mm512_load_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+48],t0r);
                              _mm512_store_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+48]);
                              jyi = _mm512_load_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+48],t1r);
                              _mm512_store_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+48]);
                              jzi = _mm512_load_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+48],t2r);
                              _mm512_store_ps(&fw.pzi[i+48],t2i);
	                 }   
	                 
	                 for(; (i+31) < n; i += 32) {
	                      cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_load_ps(&pcst[i+16]);
                              rho = _mm512_load_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+16]);
                              jxi = _mm512_load_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+16],t0r);
                              _mm512_store_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+16]);
                              jyi = _mm512_load_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+16],t1r);
                              _mm512_store_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+16]);
                              jzi = _mm512_load_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+16],t2r);
                              _mm512_store_ps(&fw.pzi[i+16],t2i);
	                 }  
	                 
	                 for(; (i+15) < n; i += 16) {
	                      cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
	                 }  
	                 
	                 for(; (i+0) < n; i += 1) {
	                     float cst               = pcst[i];
	                     float rho               = prho[i];
	                     float eai               = k*cst*rho;
	                     std::complex<float> ce = std::cexp({0.0f,eai});
	                     float jxr               = pjxr[i];
	                     float jxi               = pjxi[i];
	                     std::complex<float> cx  = {jxr,jxi};
	                     std::complex<float> jx  = cx*ce;
	                     fw.pxr[i]               = jx.real();
	                     fw.pxi[i]               = jx.imag();
	                     float jyr               = pjyr[i];
	                     float jyi               = pjyi[i];
	                     std::complex<float> cy  = {jyr,jyi};
	                     std::complex<float> jy  = cy*ce;
	                     fw.pyr[i]               = jy.real();
	                     fw.pyi[i]               = jy.imag();
	                     float jzr               = pjzr[i];
	                     float jzi               = pjzi[i];
	                     std::complex<float> cz  = {jzr,jzi};
	                     std::complex<float> jz  = ce*cz;
	                     fw.pzr[i]               = jz.real();
	                     fw.pzi[i]               = jz.imag();
	                 }                             
	       } 
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           void f2235_integrand_zmm16r4_u6x_u(const float * __restrict  pjxr,
	                                              const float * __restrict  pjxi,
	                                              const float * __restrict  pjyr,
	                                              const float * __restrict  pjyi,
	                                              const float * __restrict  pjzr,
	                                              const float * __restrict  pjzi,
	                                              const float * __restrict  prho,
	                                              const float * __restrict  pcst,
	                                              fwork_t fw, ////work arrays (caller allocated)
	                                              const float k,
	                                              const int32_t n,
	                                              const int32_t PF_DIST) {
	                                              
	                 if(__builtin_expect(n<=0,0)) { return;}
	                 register __m512 jxr,jxi,jyr,jyi,jzr,jzi;
	                 register __m512 rho,cst,ir,ii;
	                 register __m512 ear,eai,cer,cei;
	                 register __m512 vk,t0;
	                 __m512   t0r,t0i,t1r,t1i,t2r,t2i;
	                 int32_t i;
	                 
	                 vk = _mm512_set1_ps(k);
	                 ii = _mm512_set1_ps(1.0f);
	                 ir = _mm512_setzero_ps();
	                 ear= ir;
	                 t0 = _mm512_mul_ps(ii,vk); // shall be removed.
	                 for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_NTA);
#endif	                 
                              cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+16]);
                              rho = _mm512_loadu_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+16]);
                              jxi = _mm512_loadu_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+16]);
                              jyi = _mm512_loadu_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+16]);
                              jzi = _mm512_loadu_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+16],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+32]);
                              rho = _mm512_loadu_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+32]);
                              jxi = _mm512_loadu_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+32]);
                              jyi = _mm512_loadu_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+32]);
                              jzi = _mm512_loadu_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+32],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+48]);
                              rho = _mm512_loadu_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+48]);
                              jxi = _mm512_loadu_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+48]);
                              jyi = _mm512_loadu_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+48]);
                              jzi = _mm512_loadu_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+48],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+48],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+64]);
                              rho = _mm512_loadu_ps(&prho[i+64]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+64]);
                              jxi = _mm512_loadu_ps(&pjxi[i+64]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+64],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+64],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+64]);
                              jyi = _mm512_loadu_ps(&pjyi[i+64]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+64],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+64],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+64]);
                              jzi = _mm512_loadu_ps(&pjzi[i+64]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+64],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+64],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+80]);
                              rho = _mm512_loadu_ps(&prho[i+80]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+80]);
                              jxi = _mm512_loadu_ps(&pjxi[i+80]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+80],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+80],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+80]);
                              jyi = _mm512_loadu_ps(&pjyi[i+80]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+80],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+80],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+80]);
                              jzi = _mm512_loadu_ps(&pjzi[i+80]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+80],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+80],t2i);
	                 }  
	                 
	                 for(; (i+63) < n; i += 64) {
	                      cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+16]);
                              rho = _mm512_loadu_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+16]);
                              jxi = _mm512_loadu_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+16]);
                              jyi = _mm512_loadu_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+16]);
                              jzi = _mm512_loadu_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+16],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+32]);
                              rho = _mm512_loadu_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+32]);
                              jxi = _mm512_loadu_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+32]);
                              jyi = _mm512_loadu_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+32]);
                              jzi = _mm512_loadu_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+32],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+48]);
                              rho = _mm512_loadu_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+48]);
                              jxi = _mm512_loadu_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+48]);
                              jyi = _mm512_loadu_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+48]);
                              jzi = _mm512_loadu_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+48],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+48],t2i);
	                 }   
	                 
	                 for(; (i+31) < n; i += 32) {
	                       cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+16]);
                              rho = _mm512_loadu_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+16]);
                              jxi = _mm512_loadu_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+16]);
                              jyi = _mm512_loadu_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+16]);
                              jzi = _mm512_loadu_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+16],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+16],t2i);
	                 }  
	                 
	                 for(; (i+15) < n; i += 16) {
	                      cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
	                 }  
	                 
	                 for(; (i+0) < n; i += 1) {
	                     float cst               = pcst[i];
	                     float rho               = prho[i];
	                     float eai               = k*cst*rho;
	                     std::complex<float> ce = std::cexp({0.0f,eai});
	                     float jxr               = pjxr[i];
	                     float jxi               = pjxi[i];
	                     std::complex<float> cx  = {jxr,jxi};
	                     std::complex<float> jx  = cx*ce;
	                     fw.pxr[i]               = jx.real();
	                     fw.pxi[i]               = jx.imag();
	                     float jyr               = pjyr[i];
	                     float jyi               = pjyi[i];
	                     std::complex<float> cy  = {jyr,jyi};
	                     std::complex<float> jy  = cy*ce;
	                     fw.pyr[i]               = jy.real();
	                     fw.pyi[i]               = jy.imag();
	                     float jzr               = pjzr[i];
	                     float jzi               = pjzi[i];
	                     std::complex<float> cz  = {jzr,jzi};
	                     std::complex<float> jz  = ce*cz;
	                     fw.pzr[i]               = jz.real();
	                     fw.pzi[i]               = jz.imag();
	                 }                             
	       } 
	       
	       
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           void f253_integrand_zmm16r4_u6x_a(const float * __restrict __ATTR_ALIGN__(64)  phxr,
	                                             const float * __restrict __ATTR_ALIGN__(64)  phxi,
	                                             const float * __restrict __ATTR_ALIGN__(64)  phyr,
	                                             const float * __restrict __ATTR_ALIGN__(64)  phyi,
	                                             const float * __restrict __ATTR_ALIGN__(64)  phzr,
	                                             const float * __restrict __ATTR_ALIGN__(64)  phzi,
	                                             const float * __restrict __ATTR_ALIGN__(64)  pnx,
	                                             const float * __restrict __ATTR_ALIGN__(64)  pny,
	                                             const float * __restrict __ATTR_ALIGN__(64)  pnz,
	                                             float * __restrict __ATTR_ALIGN__(64)  prho,
	                                             float * __restrict __ATTR_ALIGN__(64)  pcst,
	                                             fwork_t fw, // work arrays
	                                             const float k,
	                                             const int32_t n,
	                                             const int32_t PF_DIST) {
	                                             
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m512 hxr,hxi;
	                register __m512 hyr,hyi;
	                register __m512 hzr,hzi;
	                register __m512 nx,ny,nz; 
	                register __m512 rho,cst;
	                register __m512 ear,eai;
	                register __m512 vk,t0;
	                register __m512 ir,ii;
	                         __m512 cer,cei;
	                         __m512 t0r,t0i;
	                         __m512 t1r,t1i;
	                         __m512 t2r,t2i;
	                         __m512 vxr,vxi;
	                         __m512 vyr,vyi;
	                         __m512 vzr,vzi;
	                std::complex<float> cvx,cvy,cvz;
	                int32_t i;
	                
	                vk  = _mm512_set1_ps(k);
	                cer = _mm512_setzero_ps();
	                cei = cer;
	                ii  = _mm512_set1_ps(1.0f);
	                ir  = cer;
	                ear = cer;
	                t0  = _mm512_mul_ps(ii,vk);
	                for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&phxr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&phxi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&phyr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&phyi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&phzr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&phzi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pnx[i+PF_DIST], _MM_HINT_T0);
                             _mm_prefetch((char*)&pny[i+PF_DIST], _MM_HINT_T0);
                             _mm_prefetch((char*)&pnz[i+PF_DIST], _MM_HINT_T0);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&phxr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&phxi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&phyr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&phyi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&phzr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&phzi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pnx[i+PF_DIST], _MM_HINT_T1);
                             _mm_prefetch((char*)&pny[i+PF_DIST], _MM_HINT_T1);
                             _mm_prefetch((char*)&pnz[i+PF_DIST], _MM_HINT_T1);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&phxr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&phxi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&phyr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&phyi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&phzr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&phzi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pnx[i+PF_DIST], _MM_HINT_T2);
                             _mm_prefetch((char*)&pny[i+PF_DIST], _MM_HINT_T2);
                             _mm_prefetch((char*)&pnz[i+PF_DIST], _MM_HINT_T2);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&phxr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&phxi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&phyr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&phyi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&phzr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&phzi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pnx[i+PF_DIST], _MM_HINT_NTA);
                             _mm_prefetch((char*)&pny[i+PF_DIST], _MM_HINT_NTA);
                             _mm_prefetch((char*)&pnz[i+PF_DIST], _MM_HINT_NTA);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_NTA);
#endif	                 
                             cst = _mm512_load_ps(&pcst[i+0]);
                             rho = _mm512_load_ps(&prho[i+0]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+0]);
                             hxi = _mm512_load_ps(&phxi[i+0]); 
                             hyr = _mm512_load_ps(&phyr[i+0]);
                             hyi = _mm512_load_ps(&phyi[i+0]);
                             hzr = _mm512_load_ps(&phzr[i+0]);
                             hzi = _mm512_load_ps(&phzi[i+0]);
                             nx  = _mm512_load_ps(&pnx[i+0]);
                             ny  = _mm512_load_ps(&pny[i+0]);
                             nz  = _mm512_load_ps(&pnz[i+0]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+0],t2r);
                             _mm512_store_ps(&fw.pzi[i+0],t2i);
                             cst = _mm512_load_ps(&pcst[i+16]);
                             rho = _mm512_load_ps(&prho[i+16]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+16]);
                             hxi = _mm512_load_ps(&phxi[i+16]); 
                             hyr = _mm512_load_ps(&phyr[i+16]);
                             hyi = _mm512_load_ps(&phyi[i+16]);
                             hzr = _mm512_load_ps(&phzr[i+16]);
                             hzi = _mm512_load_ps(&phzi[i+16]);
                             nx  = _mm512_load_ps(&pnx[i+16]);
                             ny  = _mm512_load_ps(&pny[i+16]);
                             nz  = _mm512_load_ps(&pnz[i+16]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+16],t2r);
                             _mm512_store_ps(&fw.pzi[i+16],t2i); 
                             cst = _mm512_load_ps(&pcst[i+32]);
                             rho = _mm512_load_ps(&prho[i+32]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+32]);
                             hxi = _mm512_load_ps(&phxi[i+32]); 
                             hyr = _mm512_load_ps(&phyr[i+32]);
                             hyi = _mm512_load_ps(&phyi[i+32]);
                             hzr = _mm512_load_ps(&phzr[i+32]);
                             hzi = _mm512_load_ps(&phzi[i+32]);
                             nx  = _mm512_load_ps(&pnx[i+32]);
                             ny  = _mm512_load_ps(&pny[i+32]);
                             nz  = _mm512_load_ps(&pnz[i+32]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+32],t0r);
                             _mm512_store_ps(&fw.pxi[i+32],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+32],t1r);
                             _mm512_store_ps(&fw.pyi[i+32],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+32],t2r);
                             _mm512_store_ps(&fw.pzi[i+32],t2i); 
                             cst = _mm512_load_ps(&pcst[i+48]);
                             rho = _mm512_load_ps(&prho[i+48]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+48]);
                             hxi = _mm512_load_ps(&phxi[i+48]); 
                             hyr = _mm512_load_ps(&phyr[i+48]);
                             hyi = _mm512_load_ps(&phyi[i+48]);
                             hzr = _mm512_load_ps(&phzr[i+48]);
                             hzi = _mm512_load_ps(&phzi[i+48]);
                             nx  = _mm512_load_ps(&pnx[i+48]);
                             ny  = _mm512_load_ps(&pny[i+48]);
                             nz  = _mm512_load_ps(&pnz[i+48]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+48],t0r);
                             _mm512_store_ps(&fw.pxi[i+48],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+48],t1r);
                             _mm512_store_ps(&fw.pyi[i+48],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+48],t2r);
                             _mm512_store_ps(&fw.pzi[i+48],t2i);
                             cst = _mm512_load_ps(&pcst[i+64]);
                             rho = _mm512_load_ps(&prho[i+64]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+64]);
                             hxi = _mm512_load_ps(&phxi[i+64]); 
                             hyr = _mm512_load_ps(&phyr[i+64]);
                             hyi = _mm512_load_ps(&phyi[i+64]);
                             hzr = _mm512_load_ps(&phzr[i+64]);
                             hzi = _mm512_load_ps(&phzi[i+64]);
                             nx  = _mm512_load_ps(&pnx[i+64]);
                             ny  = _mm512_load_ps(&pny[i+64]);
                             nz  = _mm512_load_ps(&pnz[i+64]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+64],t0r);
                             _mm512_store_ps(&fw.pxi[i+64],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+64],t1r);
                             _mm512_store_ps(&fw.pyi[i+64],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+64],t2r);
                             _mm512_store_ps(&fw.pzi[i+64],t2i);
                             cst = _mm512_load_ps(&pcst[i+80]);
                             rho = _mm512_load_ps(&prho[i+80]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+80]);
                             hxi = _mm512_load_ps(&phxi[i+80]); 
                             hyr = _mm512_load_ps(&phyr[i+80]);
                             hyi = _mm512_load_ps(&phyi[i+80]);
                             hzr = _mm512_load_ps(&phzr[i+80]);
                             hzi = _mm512_load_ps(&phzi[i+80]);
                             nx  = _mm512_load_ps(&pnx[i+80]);
                             ny  = _mm512_load_ps(&pny[i+80]);
                             nz  = _mm512_load_ps(&pnz[i+80]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+80],t0r);
                             _mm512_store_ps(&fw.pxi[i+80],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+80],t1r);
                             _mm512_store_ps(&fw.pyi[i+80],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+80],t2r);
                             _mm512_store_ps(&fw.pzi[i+80],t2i);	                
	                }   
	                
	                for(; (i+63) < n; i += 64) {
	                     cst = _mm512_load_ps(&pcst[i+0]);
                             rho = _mm512_load_ps(&prho[i+0]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+0]);
                             hxi = _mm512_load_ps(&phxi[i+0]); 
                             hyr = _mm512_load_ps(&phyr[i+0]);
                             hyi = _mm512_load_ps(&phyi[i+0]);
                             hzr = _mm512_load_ps(&phzr[i+0]);
                             hzi = _mm512_load_ps(&phzi[i+0]);
                             nx  = _mm512_load_ps(&pnx[i+0]);
                             ny  = _mm512_load_ps(&pny[i+0]);
                             nz  = _mm512_load_ps(&pnz[i+0]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+0],t2r);
                             _mm512_store_ps(&fw.pzi[i+0],t2i);
                             cst = _mm512_load_ps(&pcst[i+16]);
                             rho = _mm512_load_ps(&prho[i+16]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+16]);
                             hxi = _mm512_load_ps(&phxi[i+16]); 
                             hyr = _mm512_load_ps(&phyr[i+16]);
                             hyi = _mm512_load_ps(&phyi[i+16]);
                             hzr = _mm512_load_ps(&phzr[i+16]);
                             hzi = _mm512_load_ps(&phzi[i+16]);
                             nx  = _mm512_load_ps(&pnx[i+16]);
                             ny  = _mm512_load_ps(&pny[i+16]);
                             nz  = _mm512_load_ps(&pnz[i+16]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+16],t2r);
                             _mm512_store_ps(&fw.pzi[i+16],t2i); 
                             cst = _mm512_load_ps(&pcst[i+32]);
                             rho = _mm512_load_ps(&prho[i+32]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+32]);
                             hxi = _mm512_load_ps(&phxi[i+32]); 
                             hyr = _mm512_load_ps(&phyr[i+32]);
                             hyi = _mm512_load_ps(&phyi[i+32]);
                             hzr = _mm512_load_ps(&phzr[i+32]);
                             hzi = _mm512_load_ps(&phzi[i+32]);
                             nx  = _mm512_load_ps(&pnx[i+32]);
                             ny  = _mm512_load_ps(&pny[i+32]);
                             nz  = _mm512_load_ps(&pnz[i+32]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+32],t0r);
                             _mm512_store_ps(&fw.pxi[i+32],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+32],t1r);
                             _mm512_store_ps(&fw.pyi[i+32],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+32],t2r);
                             _mm512_store_ps(&fw.pzi[i+32],t2i); 
                             cst = _mm512_load_ps(&pcst[i+48]);
                             rho = _mm512_load_ps(&prho[i+48]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+48]);
                             hxi = _mm512_load_ps(&phxi[i+48]); 
                             hyr = _mm512_load_ps(&phyr[i+48]);
                             hyi = _mm512_load_ps(&phyi[i+48]);
                             hzr = _mm512_load_ps(&phzr[i+48]);
                             hzi = _mm512_load_ps(&phzi[i+48]);
                             nx  = _mm512_load_ps(&pnx[i+48]);
                             ny  = _mm512_load_ps(&pny[i+48]);
                             nz  = _mm512_load_ps(&pnz[i+48]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+48],t0r);
                             _mm512_store_ps(&fw.pxi[i+48],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+48],t1r);
                             _mm512_store_ps(&fw.pyi[i+48],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+48],t2r);
                             _mm512_store_ps(&fw.pzi[i+48],t2i);
	                }  
	                
	                for(; (i+31) < n; i += 32) {
	                     cst = _mm512_load_ps(&pcst[i+0]);
                             rho = _mm512_load_ps(&prho[i+0]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+0]);
                             hxi = _mm512_load_ps(&phxi[i+0]); 
                             hyr = _mm512_load_ps(&phyr[i+0]);
                             hyi = _mm512_load_ps(&phyi[i+0]);
                             hzr = _mm512_load_ps(&phzr[i+0]);
                             hzi = _mm512_load_ps(&phzi[i+0]);
                             nx  = _mm512_load_ps(&pnx[i+0]);
                             ny  = _mm512_load_ps(&pny[i+0]);
                             nz  = _mm512_load_ps(&pnz[i+0]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+0],t2r);
                             _mm512_store_ps(&fw.pzi[i+0],t2i);
                             cst = _mm512_load_ps(&pcst[i+16]);
                             rho = _mm512_load_ps(&prho[i+16]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+16]);
                             hxi = _mm512_load_ps(&phxi[i+16]); 
                             hyr = _mm512_load_ps(&phyr[i+16]);
                             hyi = _mm512_load_ps(&phyi[i+16]);
                             hzr = _mm512_load_ps(&phzr[i+16]);
                             hzi = _mm512_load_ps(&phzi[i+16]);
                             nx  = _mm512_load_ps(&pnx[i+16]);
                             ny  = _mm512_load_ps(&pny[i+16]);
                             nz  = _mm512_load_ps(&pnz[i+16]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+16],t2r);
                             _mm512_store_ps(&fw.pzi[i+16],t2i); 
	                }
	                
	                for(; (i+15) < n; i += 16) {
	                     cst = _mm512_load_ps(&pcst[i+0]);
                             rho = _mm512_load_ps(&prho[i+0]);
                             eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                             cexp_zmm16r4(ear,eai,&cer,&cei); 
                             hxr = _mm512_load_ps(&phxr[i+0]);
                             hxi = _mm512_load_ps(&phxi[i+0]); 
                             hyr = _mm512_load_ps(&phyr[i+0]);
                             hyi = _mm512_load_ps(&phyi[i+0]);
                             hzr = _mm512_load_ps(&phzr[i+0]);
                             hzi = _mm512_load_ps(&phzi[i+0]);
                             nx  = _mm512_load_ps(&pnx[i+0]);
                             ny  = _mm512_load_ps(&pny[i+0]);
                             nz  = _mm512_load_ps(&pnz[i+0]);
                             scrosscv_zmm16c4(hxr,hxi,hyr,
                                              hyi,hzr,hzi,
                                              nx,ny,nz,
                                              &vxr,&vxi,&vyr,
                                              &vyi,&vzr,&vzi);
                             cmul_zmm16r4(vxr,vxi,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             cmul_zmm16r4(vyr,vyi,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             cmul_zmm16r4(vzr,vzi,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.pzr[i+0],t2r);
                             _mm512_store_ps(&fw.pzi[i+0],t2i);
	                }
	                
	                
	                for(; (i+0) < n; i += 1) {
	                     float cst               = pcst[i];
	                     float rho               = prho[i];
	                     float eai               = k*cst*rho;
	                     std::complex<float> ce  = std::cexp({0.0f,eai});
	                     float hxr               = phxr[i];
	                     float hxi               = phxi[i];
	                     float hyr               = phyr[i];
	                     float hyi               = phyi[i];
	                     float hzr               = phzr[i];
	                     float hzi               = phzi[i];
	                     float nx                = pnx[i];
	                     float ny                = pny[i];
	                     float nz                = pnz[i];
	                     scrosscv_rc4(hxr,hxi,hyr,
	                                  hyi,hzr,hzi,
	                                  nx,ny,nz,
	                                  cvx,cvy,cvz);
	                     std::complex<float> t0  = ce*cvx;
	                     fw.pxr[i]               = t0.real();
	                     fw.pxi[i]               = t0.imag();
	                     std::complex<float> t1  = ce*cvy;
	                     fw.pyr[i]               = t1.real();
	                     fw.pyi[i]               = t1.imag();
	                     std::complex<float> t2  = ce*cvz;
	                     fw.pzr[i]               = t2.real();
	                     fw.pzi[i]               = t2.imag();
	                }
	                                                     
	     }
	       
	       
	        
                  
                  
       
       } // radiolocation
       
       
} // gms




#endif /*__GMS_ANTENNA_FEEDER_INTEGRANDS_ZMM16R4*/
