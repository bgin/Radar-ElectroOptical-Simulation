
#ifndef __GMS_COMPLEX_XMM4R4_HPP__
#define __GMS_COMPLEX_XMM4R4_HPP__ 22120231700


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

    const unsigned int GMS_COMPLEX_XMM4R4_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_XMM4R4_MINOR = 0U;
    const unsigned int GMS_COMPLEX_XMM4R4_MICRO = 0U;
    const unsigned int GMS_COMPLEX_XMM4R4_FULLVER =
      1000U*GMS_COMPLEX_XMM4R4_MAJOR+
      100U*GMS_COMPLEX_XMM4R4_MINOR+
      10U*GMS_COMPLEX_XMM4R4_MICRO;
    const char * const GMS_COMPLEX_XMM4R4_CREATION_DATE = "22-10-2023 17:00  +00200 (SUN 22 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_XMM4R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_XMM4R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_XMM4R4_DESCRIPTION   = "SSE optimized complex number implementation."

}

#include <cstdint>
#include <eimmintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.hpp"

namespace  gms {


       namespace math {
       
       
                   struct __ATTR_ALIGN__(16) xmm4c4_t {
                   
                          __m128 re;
                          __m128 im;
                   };
                   
                   

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_ps(&xre[0]);
                        xmm1  = _mm_loadu_ps(&yre[0]);
                        _mm_storeu_ps(&zre[0], _mm_add_ps(xmm0,xmm1));
                        xmm2  = _mm_loadu_ps(&xim[0]);
                        xmm3  = _mm_loadu_ps(&yim[0]);
                        _mm_storeu_ps(&zim[0], _mm_add_ps(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       const float * __restrict __ATTR_ALIGN__(16) yre,
                                       const float * __restrict __ATTR_ALIGN__(16) yim,
                                       float *       __restrict __ATTR_ALIGN__(16) zre,
                                       float *       __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_ps(&xre[0]);
                        xmm1  = _mm_load_ps(&yre[0]);
                        _mm_store_ps(&zre[0], _mm_add_ps(xmm0,xmm1));
                        xmm2  = _mm_load_ps(&xim[0]);
                        xmm3  = _mm_load_ps(&yim[0]);
                        _mm_store_ps(&zim[0], _mm_add_ps(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 yre,
                                     const __m128 yim,
                                     __m128 * __restrict zre,
                                     __m128 * __restrict zim) {
                     
                        register __m128 xmm0,xmm1;
                        xmm0  = _mm_add_ps(xre,yre);
                        *zre  = xmm0;
                        xmm1  = _mm_add_ps(xim,yim);
                        *zim  = xmm1;
                }
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cadd_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 yre,
                                          const __m128 yim) {
                                     
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0   = _mm_add_ps(xre,yre);
                        cv.re  = xmm0;
                        xmm1   = _mm_add_ps(xim,yim);
                        cv.im  = xmm1;  
                        return (cv);            
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cadd_xmm4c4(const xmm4c4_t x,
                                          const xmm4c4_t y) {
                                     
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0   = _mm_add_ps(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1   = _mm_add_ps(x.im,y.im);
                        cv.im  = xmm1;  
                        return (cv);            
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 s,
                                     __m128 * __restrict     zre,
                                     __m128 * __restrict     zim) {

                        *zre = _mm_add_ps(xre,s);
                        *zim = _mm_setzero_ps();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           xmm4c4_t cadd_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 s) {
                      
                      xmm4c4_t cv;
                      cv.re =  _mm_add_ps(xre,s);
                      cv.im =  _mm_setzero_ps();
                      return (cv);                       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           xmm4c4_t cadd_xmm4c4(const xmm4c4_t x,
                                          const __m128 s) {
                      
                      xmm4c4_t cv;
                      cv.re =  _mm_add_ps(x.re,s);
                      cv.im =  _mm_setzero_ps();
                      return (cv);                       
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_xmm4c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_ps(&xre[0]);
                        xmm1  = _mm_loadu_ps(&xim[0]);
                        xmm2  = _mm_loadu_ps(&zre[0]);
                        xmm3  = _mm_loadu_ps(&zim[0])
                        _mm_storeu_ps(&zre[0], _mm_add_ps(xmm2,xmm0));
                        _mm_storeu_ps(&zim[0], _mm_add_ps(xmm3,xmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_xmm4c4_aip(const float * __restrict __ATTR_ALIGN__(16) xre,
                                         const float * __restrict __ATTR_ALIGN__(16) xim,
                                         float *       __restrict __ATTR_ALIGN__(16) zre,
                                         float *       __restrict __ATTR_ALIGN__(16) zim) {
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_ps(&xre[0]);
                        xmm1  = _mm_load_ps(&xim[0]);
                        xmm2  = _mm_load_ps(&zre[0]);
                        xmm3  = _mm_load_ps(&zim[0])
                        _mm_store_ps(&zre[0], _mm_add_ps(xmm2,xmm0));
                        _mm_store_ps(&zim[0], _mm_add_ps(xmm3,xmm1));
              }


                ////////////////////////////////////////////////////////////////////

                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_ps(&xre[0]);
                        xmm1  = _mm_loadu_ps(&yre[0]);
                        _mm_storeu_ps(&zre[0], _mm_sub_ps(xmm0,xmm1));
                        xmm2  = _mm_loadu_ps(&xim[0]);
                        xmm3  = _mm_loadu_ps(&yim[0]);
                        _mm_storeu_ps(&zim[0], _mm_sub_ps(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_xmm4c4_a( const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       const float * __restrict __ATTR_ALIGN__(16) yre,
                                       const float * __restrict __ATTR_ALIGN__(16) yim,
                                       float *       __restrict __ATTR_ALIGN__(16) zre,
                                       float *       __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_ps(&xre[0]);
                        xmm1  = _mm_load_ps(&yre[0]);
                        _mm_store_ps(&zre[0], _mm_sub_ps(xmm0,xmm1));
                        xmm2  = _mm_load_ps(&xim[0]);
                        xmm3  = _mm_load_ps(&yim[0]);
                        _mm_store_ps(&zim[0], _mm_sub_ps(xmm2,xmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 yre,
                                     const __m128 yim,
                                     __m128 * __restrict  zre,
                                     __m128 * __restrict  zim) {
                     
                        register __m128 xmm0,xmm1;
                        xmm0  = _mm_sub_ps(xre,yre);
                        *zre  = xmm0;
                        xmm1  = _mm_sub_ps(xim,yim);
                        *zim  = xmm1;
                }
                
                
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t csub_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 yre,
                                          const __m128 yim) {
                                    
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0  = _mm_sub_ps(xre,yre);
                        cv.re  = xmm0;
                        xmm1  = _mm_sub_ps(xim,yim);
                        cv.im  = xmm1;
                        return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t csub_xmm4c4(const xmm4c4_t x,
                                          const xmm4c4_t y) {
                                    
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0  = _mm_sub_ps(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1  = _mm_sub_ps(x.im,y.im);
                        cv.im  = xmm1;
                        return (cv);
                }
                
                


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 s,
                                     __m128 * __restrict     zre,
                                     __m128 * __restrict     zim) {

                        *zre = _mm_sub_ps(xre,s);
                        *zim = _mm_setzero_ps();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t csub_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 s) {
                                    
                        xmm4c4_t cv;
                        cv.re = _mm_sub_ps(xre,s);
                        cv.im = _mm_setzero_ps();
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t csub_xmm4c4(const xmm4c4_t x,
                                          const __m128 s) {
                                    
                        xmm4c4_t cv;
                        cv.re = _mm_sub_ps(x.re,s);
                        cv.im = _mm_setzero_ps();
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_xmm4c4_uip( const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_ps(&xre[0]);
                        xmm1  = _mm_loadu_ps(&xim[0]);
                        xmm2  = _mm_loadu_ps(&zre[0]);
                        xmm3  = _mm_loadu_ps(&zim[0])
                        _mm_storeu_ps(&zre[0], _mm_sub_ps(xmm2,xmm0));
                        _mm_storeu_ps(&zim[0], _mm_sub_ps(xmm3,xmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_xmm4c4_aip(const float * __restrict __ATTR_ALIGN__(16) xre,
                                         const float * __restrict __ATTR_ALIGN__(16) xim,
                                         float *       __restrict __ATTR_ALIGN__(16) zre,
                                         float *       __restrict __ATTR_ALIGN__(16) zim) {
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_ps(&xre[0]);
                        xmm1  = _mm_load_ps(&xim[0]);
                        xmm2  = _mm_load_ps(&zre[0]);
                        xmm3  = _mm_load_ps(&zim[0])
                        _mm_store_ps(&zre[0], _mm_sub_ps(xmm2,xmm0));
                        _mm_store_ps(&zim[0], _mm_sub_ps(xmm3,xmm1));
              }


               ////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                           register __m128 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_loadu_ps(&xre[0]);
                           xmm1  = _mm_loadu_ps(&yre[0]);
                           xmm2  = _mm_loadu_ps(&xim[0]);
                           xmm3  = _mm_loadu_ps(&yim[0]);
                           xmm4  = _mm_sub_ps(_mm_mul_ps(xmm0,xmm1),
                                                   _mm_mul_ps(xmm2,xmm3));
                           _mm_storeu_ps(&zre[0], xmm4);
                           xmm5  = _mm_mul_ps(_mm_mul_ps(xmm2,xmm1),
                                                   _mm_mul_ps(xmm0,xmm3));
                           _mm_storeu_ps(&zim[0], xmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       const float * __restrict __ATTR_ALIGN__(16) yre,
                                       const float * __restrict __ATTR_ALIGN__(16) yim,
                                       float *       __restrict __ATTR_ALIGN__(16) zre,
                                       float *       __restrict __ATTR_ALIGN__(16) zim) {

                           register __m128 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_load_ps(&xre[0]);
                           xmm1  = _mm_load_ps(&yre[0]);
                           xmm2  = _mm_load_ps(&xim[0]);
                           xmm3  = _mm_load_ps(&yim[0]);
                           xmm4  = _mm_sub_ps(_mm_mul_ps(xmm0,xmm1),
                                                                        _mm_mul_ps(xmm2,xmm3));
                           _mm_store_ps(&zre[0], xmm4);
                           xmm5  = _mm_mul_ps(_mm_mul_ps(xmm2,xmm1),
                                                                        _mm_mul_ps(xmm0,xmm3));
                           _mm_store_ps(&zim[0], xmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 yre,
                                     const __m128 yim,
                                     __m128 * __restrict     zre,
                                     __m128 * __restrict     zim) {

                         register __m128 xmm0,xmm1;
                         xmm0 = _mm_sub_ps(_mm_mul_ps(xre,yre),
                                              _mm_mul_ps(xim,yim));
                         *zre  = xmm0;
                         xmm1 = _mm_mul_ps(_mm_mul_ps(xim,yre),
                                              _mm_mul_ps(xre,yim));
                         *zim  = xmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cmul_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 yre,
                                          const __m128 yim) {
                                     
                         xmm4c4_t cv
                         register __m128 xmm0,xmm1;
                         xmm0 = _mm_sub_ps(_mm_mul_ps(xre,yre),
                                              _mm_mul_ps(xim,yim));
                         cv.re  = xmm0;
                         xmm1 = _mm_mul_ps(_mm_mul_ps(xim,yre),
                                              _mm_mul_ps(xre,yim));
                         cv.im  = xmm1;
                         return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cmul_xmm4c4(const xmm4c4_t x,
                                          const xmm4c4_t y) {
                                     
                         xmm4c4_t cv
                         register __m128 xmm0,xmm1;
                         xmm0 = _mm_sub_ps(_mm_mul_ps(x.re,y.re),
                                              _mm_mul_ps(x.im,y.im));
                         cv.re  = xmm0;
                         xmm1 = _mm_mul_ps(_mm_mul_ps(x.im,y.re),
                                              _mm_mul_ps(x.re,y.im));
                         cv.im  = xmm1;
                         return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 s,
                                     __m128 * __restrict   zre,
                                     __m128 * __restrict   zim) {

                        *zre = _mm_mul_ps(xre,s);
                        *zim = _mm_mul_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cmul_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 s) {
                                     
                        xmm4c4_t cv;
                        cv.re = _mm_mul_ps(xre,s);
                        cv.im = _mm_mul_ps(xim,s);
                        return (cv);
               }
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cmul_xmm4c4(const xmm4c4_t x,
                                          const __m128 s) {
                                     
                        xmm4c4_t cv;
                        cv.re = _mm_mul_ps(x.re,s);
                        cv.im = _mm_mul_ps(x.im,s);
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_xmm4c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {

                           register __m128 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_loadu_ps(&xre[0]);
                           xmm1  = _mm_loadu_ps(&zre[0]);
                           xmm2  = _mm_loadu_ps(&xim[0]);
                           xmm3  = _mm_loadu_ps(&zim[0]);
                           xmm4  = _mm_sub_ps(_mm_mul_ps(xmm0,xmm1),
                                                 _mm_mul_ps(xmm2,xmm3));
                           _mm_storeu_ps(&zre[0], xmm4);
                           xmm5  = _mm_mul_ps(_mm_mul_ps(xmm2,xmm1),
                                                 _mm_mul_ps(xmm0,xmm3));
                           _mm_storeu_ps(&zim[0], xmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_xmm4c4_aip(const float * __restrict __ATTR_ALIGN__(16) xre,
                                         const float * __restrict __ATTR_ALIGN__(16) xim,
                                         float *       __restrict __ATTR_ALIGN__(16) zre,
                                         float *       __restrict __ATTR_ALIGN__(16) zim) {

                           register __m128 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
                           xmm0  = _mm_load_ps(&xre[0]);
                           xmm1  = _mm_load_ps(&zre[0]);
                           xmm2  = _mm_load_ps(&xim[0]);
                           xmm3  = _mm_load_ps(&zim[0]);
                           xmm4  = _mm_sub_ps(_mm_mul_ps(xmm0,xmm1),
                                                 _mm_mul_ps(xmm2,xmm3));
                           _mm_store_ps(&zre[0], xmm4);
                           xmm5  = _mm_mul_ps(_mm_mul_ps(xmm2,xmm1),
                                                 _mm_mul_ps(xmm0,xmm3));
                           _mm_store_ps(&zim[0], xmm5);
               }

                 ////////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3; 
                        register __m128 xmm4,xmm5,xmm6;
                        xmm0  = _mm_loadu_ps(&xre[0]); //a
                        xmm1  = _mm_loadu_ps(&yim[0]); //d
                        xmm2  = _mm_loadu_ps(&xim[0]); //b
                        xmm3  = _mm_loadu_ps(&yre[0]); //c
                        xmm4  = _mm_fmadd_ps(xmm0,xmm3,
                                                _mm_mul_ps(xmm2,xmm1));
                        xmm5  = _mm_fmsub_ps(xmm2,xmm3,
                                                _mm_mul_ps(xmm0,xmm1));
                        xmm6  = _mm_fmadd_ps(xmm3,xmm3),
                                                _mm_mul_ps(xmm1,xmm1));
                        _mm_storeu_ps(&zre[0], _mm_div_ps(xmm4,xmm6));
                        _mm_storeu_ps(&zim[0], _mm_div_ps(xmm5,xmm6));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4_a(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3; 
                        register __m128 xmm4,xmm5,xmm6;
                        xmm0  = _mm_load_ps(&xre[0]); //a
                        xmm1  = _mm_load_ps(&yim[0]); //d
                        xmm2  = _mm_load_ps(&xim[0]); //b
                        xmm3  = _mm_load_ps(&yre[0]); //c
                        xmm4  = _mm_fmadd_ps(xmm0,xmm3,
                                                _mm_mul_ps(xmm2,xmm1));
                        xmm5  = _mm_fmsub_ps(xmm2,xmm3,
                                                _mm_mul_ps(xmm0,xmm1));
                        xmm6  = _mm_fmadd_ps(xmm3,xmm3,
                                                _mm_mul_ps(xmm1,xmm1));
                        _mm_store_ps(&zre[0], _mm_div_ps(xmm4,xmm6));
                        _mm_store_ps(&zim[0], _mm_div_ps(xmm5,xmm6));
                }
              
                                         
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 yre,
                                     const __m128 yim,
                                     __m128 * __restrict zre,
                                     __m128 * __restrict zim) {

                      register __m128 xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_ps(xre,yre,
                                           _mm_mul_ps(xim,yim));
                      xmm1 = _mm_fmsub_ps(xim,yre,
                                           _mm_mul_ps(xre,yim));
                      xmm2 = _mm_fmadd_ps(xmm3,xmm3,
                                           _mm_mul_ps(xmm1,xmm1));
                      *zre  = _mm_div_ps(xmm0,xmm2);
                      *zim  = _mm_div_ps(xmm1,xmm2);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 yre,
                                          const __m128 yim) {
                                     
                      xmm4c4_t
                      register __m128 xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_ps(xre,yre,
                                           _mm_mul_ps(xim,yim));
                      xmm1 = _mm_fmsub_ps(x.im,y.re,
                                           _mm_mul_ps(xre,yim));
                      xmm2 = _mm_fmadd_ps(xmm3,xmm3,
                                           _mm_mul_ps(xmm1,xmm1));
                      cv.re  = _mm_div_ps(xmm0,xmm2);
                      cv.im  = _mm_div_ps(xmm1,xmm2);
                      return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_xmm4c4(const xmm4c4_t x,
                                          const xmm4c4_t y) {
                                     
                      xmm4c4_t
                      register __m128 xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_ps(x.re,y.re,
                                           _mm_mul_ps(x.im,y.im));
                      xmm1 = _mm_fmsub_ps(x.im,y.re,
                                           _mm_mul_ps(x.re,y.im));
                      xmm2 = _mm_fmadd_ps(xmm3,xmm3,
                                           _mm_mul_ps(xmm1,xmm1));
                      cv.re  = _mm_div_ps(xmm0,xmm2);
                      cv.im  = _mm_div_ps(xmm1,xmm2);
                      return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     const __m128 s,
                                     __m128 * __restrict zre,
                                     __m128 * __restrict zim) {

                        *zre = _mm_div_ps(xre,s);
                        *zim = _mm_div_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_xmm4c4(const __m128 xre,
                                          const __m128 xim,
                                          const __m128 s) {
                                     
                         xmm4c4_t cv;
                         cv.re = _mm_div_ps(xre,s);
                         cv.im = _mm_div_ps(xim,s);
                         return (cv);
               }
               
               
                 
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_xmm4c4(const xmm4c4_t x,
                                          const __m128 s) {
                                     
                         xmm4c4_t cv;
                         cv.re = _mm_div_ps(x.re,s);
                         cv.im = _mm_div_ps(x.im,s);
                         return (cv);
               }
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4_s(const __m128 s,
                                       const __m128 xre,
                                       const __m128 xim,
                                       __m128 * __restrict zre,
                                       __m128 * __restrict zim) {
                        
                        register __m128 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm_setzero_ps();
                        cdiv_xmm4c4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_xmm4c4_s(const __m128 s,
                                            const __m128 xre,
                                            const __m128 xim) {
                                       
                        xmm4c4_t cv;
                        register __m128 t0r,t0i;
                        t0r = s;
                        t0i = _mm_setzero_ps();
                        cdiv_xmm4c4(t0r,t0i,xre,xim,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_xmm4c4_s(const __m128 s,
                                            const xmm4c4_t x) {
                                       
                        xmm4c4_t cv;
                        register __m128 t0r,t0i;
                        t0r = s;
                        t0i = _mm_setzero_ps();
                        cdiv_xmm4c4(t0r,t0i,x.re,x.im,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_xmm4c4_s(const __m128 s,
                                             const __m128 xre,
                                             const __m128 xim,
                                             __m128 * __restrict zre,
                                             __m128 * __restrict zim) {
                                             
                        register __m128 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm_setzero_ps(); 
                        cdiv_smith_xmm4c4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                   
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_smith_xmm4c4_s(const __m128 s,
                                                  const __m128 xre,
                                                  const __m128 xim) {
                                             
                        xmm4c4_t cv;                    
                        register __m128 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm_setzero_ps(); 
                        cdiv_smith_xmm4c4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_smith_xmm4c4_s(const __m128 s,
                                                  const xmm4c4_t x) {
                                             
                        xmm4c4_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm_setzero_ps(); 
                        cv = cdiv_smith_xmm4c4(t0,x);
                        return (cv);                 
                 }
                 
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3; 
                        register __m128 xmm4,xmm5,xmm6;
                        xmm0  = _mm_loadu_ps(&xre[0]); //a
                        xmm1  = _mm_loadu_ps(&zim[0]); //d
                        xmm2  = _mm_loadu_ps(&xim[0]); //b
                        xmm3  = _mm_loadu_ps(&zre[0]); //c
                        xmm4  = _mm_fmadd_ps(xmm0,xmm3,
                                                _mm_mul_ps(xmm2,xmm1));
                        xmm5  = _mm_fmsub_ps(xmm2,xmm3,
                                                _mm_mul_ps(xmm0,xmm1));
                        xmm6  = _mm_fmadd_ps(xmm3,xmm3,
                                                _mm_mul_ps(xmm1,xmm1));
                        _mm_storeu_ps(&zre[0], _mm_div_ps(xmm4,xmm6));
                        _mm_storeu_ps(&zim[0], _mm_div_ps(xmm5,xmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_xmm4c4_aip(const float * __restrict __ATTR_ALIGN__(16) xre,
                                         const float * __restrict __ATTR_ALIGN__(16) xim,
                                         float *       __restrict __ATTR_ALIGN__(16) zre,
                                         float *       __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3; 
                        register __m128 xmm4,xmm5,xmm6;
                        xmm0  = _mm_load_ps(&xre[0]); //a
                        xmm1  = _mm_load_ps(&zim[0]); //d
                        xmm2  = _mm_load_ps(&xim[0]); //b
                        xmm3  = _mm_load_ps(&zre[0]); //c
                        xmm4  = _mm_fmadd_ps(xmm0,xmm3,
                                                _mm_mul_ps(xmm2,xmm1));
                        xmm5  = _mm_fmsub_ps(xmm2,xmm3,
                                                _mm_mul_ps(xmm0,xmm1));
                        xmm6  = _mm_fmadd_ps(xmm3,xmm3,
                                                _mm_mul_ps(xmm1,xmm1));
                        _mm_store_ps(&zre[0], _mm_div_ps(xmm4,xmm6));
                        _mm_store_ps(&zim[0], _mm_div_ps(xmm5,xmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_xmm4c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float *       __restrict zre,
                                             float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 r,den;
                        __mmask8 m = 0x0;
                        xmm0 = _mm_loadu_ps(&yre[0]); // c
                        xmm1 = _mm_loadu_ps(&yim[0]); // d
                        xmm2 = _mm_loadu_ps(&xre[0]); // a
                        xmm3 = _mm_loadu_ps(&xim[0]); // b
                        m    = _mm_cmp_ps_mask(_mm_abs_ps(xmm0),
                                                  _mm_abs_ps(xmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_ps(m,_mm_div_ps(xmm0,xmm1),
                                                      _mm_div_ps(xmm1,xmm0)); // r
                        den  = _mm_mask_blend_ps(m,_mm_fmadd_ps(r,xmm0,xmm1),
                                                      _mm_fmadd_ps(r,xmm1,xmm0));
                        _mm_storeu_ps(&zre[0], _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmadd_ps(xmm2,r,xmm3),den),
                                                _mm_div_ps(_mm_fmadd_ps(xmm3,r,xmm2),den)));
                        _mm_storeu_ps(&zim[0], _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmsub_ps(xmm3,r,xmm2),den),
                                                _mm_div_ps(_mm_sub_ps(xmm3,_mm_mul_ps(xmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_xmm4c4_a(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float *       __restrict zre,
                                             float *       __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 r,den;
                        __mmask8 m = 0x0;
                        xmm0 = _mm_load_ps(&yre[0]); // c
                        xmm1 = _mm_load_ps(&yim[0]); // d
                        xmm2 = _mm_load_ps(&xre[0]); // a
                        xmm3 = _mm_load_ps(&xim[0]); // b
                        m    = _mm_cmp_ps_mask(_mm_abs_ps(xmm0),
                                                  _mm_abs_ps(xmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_ps(m,_mm_div_ps(xmm0,xmm1),
                                                      _mm_div_ps(xmm1,xmm0)); // r
                        den  = _mm_mask_blend_ps(m,_mm_fmadd_ps(r,xmm0,xmm1),
                                                      _mm_fmadd_ps(r,xmm1,xmm0));
                        _mm_storeu_ps(&zre[0], _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmadd_ps(xmm2,r,xmm3),den),
                                                _mm_div_ps(_mm_fmadd_ps(xmm3,r,xmm2),den)));
                        _mm_storeu_ps(&zim[0], _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmsub_ps(xmm3,r,xmm2),den),
                                                _mm_div_ps(_mm_sub_ps(xmm3,_mm_mul_ps(xmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_xmm4c4(const __m128 xre,
                                           const __m128 xim,
                                           const __m128 yre,
                                           const __m128 yim,
                                           __m128 * __restrict zre,
                                           __m128 * __restrict zim) {

                        register __m128 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_ps_mask(_mm_abs_ps(yre),
                                                  _mm_abs_ps(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_ps(m,_mm_div_ps(yre,yim),
                                                      _mm_div_ps(yim,yre)); // r
                        den  = _mm_mask_blend_ps(m,_mm_fmadd_ps(r,yre,yim),
                                                      _mm_fmadd_ps(r,yim,yre));
                        *zre  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmadd_ps(xre,r,xim),den),
                                                _mm_div_ps(_mm_fmadd_ps(xim,r,xre),den));
                        *zim  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmsub_ps(xim,r,xre),den),
                                                _mm_div_ps(_mm_sub_ps(xim,_mm_mul_ps(xre,r)),den)));
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_smith_xmm4c4(const __m128 xre,
                                                const __m128 xim,
                                                const __m128 yre,
                                                const __m128 yim) {
                                           
                        xmm4c4_t cv
                        register __m128 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_ps_mask(_mm_abs_ps(yre),
                                                  _mm_abs_ps(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_ps(m,_mm_div_ps(yre,yim),
                                                      _mm_div_ps(yim,yre)); // r
                        den  = _mm_mask_blend_ps(m,_mm_fmadd_ps(r,yre,yim),
                                                      _mm_fmadd_ps(r,yim,yre));
                        cv.re  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmadd_ps(xre,r,xim),den),
                                                _mm_div_ps(_mm_fmadd_ps(xim,r,xre),den));
                        cv.im  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmsub_ps(xim,r,xre),den),
                                                _mm_div_ps(_mm_sub_ps(xim,_mm_mul_ps(xre,r)),den)));
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cdiv_smith_xmm4c4(const xmm4c4_t x,
                                                const xmm4c4_t y) {
                                           
                        xmm4c4_t cv
                        register __m128 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_ps_mask(_mm_abs_ps(y.re),
                                                  _mm_abs_ps(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_ps(m,_mm_div_ps(y.re,y.im),
                                                      _mm_div_ps(y.im,y.re)); // r
                        den  = _mm_mask_blend_ps(m,_mm_fmadd_ps(r,y.re,y.im),
                                                      _mm_fmadd_ps(r,y.im,y.re));
                        cv.re  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmadd_ps(x.re,r,x.im),den),
                                                _mm_div_ps(_mm_fmadd_ps(x.im,r,x.re),den));
                        cv.im  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmsub_ps(x.im,r,x.re),den),
                                                _mm_div_ps(_mm_sub_ps(x.im,_mm_mul_ps(x.re,r)),den)));
                        return (cv);
               }





                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_xmm4c4_u(const float * __restrict re,
                                       const float * __restrict im,
                                       float * __restrict  cabs) {

                        register __m128 xmm0,xmm1,xmm2,xmm3,xmm4;
                        xmm0  = _mm_loadu_ps(&re[0]);
                        xmm1  = _mm_mul_ps(xmm0,xmm0);
                        xmm2  = _mm_loadu_ps(&im[0]);
                        xmm3  = _mm_mul_ps(xmm2,xmm2);
                        xmm4  = _mm_sqrt_ps(_mm_add_ps(xmm1,xmm3));
                        _mm_storeu_ps(&cabs[0],xmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) re,
                                       const float * __restrict __ATTR_ALIGN__(16) im,
                                       float * __restrict  __ATTR_ALIGN__(16) cabs) {

                        register __m128 xmm0,xmm1,xmm2,xmm3,xmm4;
                        xmm0  = _mm_load_ps(&re[0]);
                        xmm1  = _mm_mul_ps(xmm0,xmm0);
                        xmm2  = _mm_load_ps(&im[0]);
                        xmm3  = _mm_mul_ps(xmm2,xmm2);
                        xmm4  = _mm_sqrt_ps(_mm_add_ps(xmm1,xmm3));
                        _mm_store_ps(&cabs[0],xmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m128 cabs_xmm4c4(const __m128 re,
                                       const __m128 im) {

                        register __m128 xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_ps(re,re);
                        xmm1 = _mm_mul_ps(im,im);
                        cabs = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        return (cabs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m128 cabs_xmm4c4(const xmm4c4_t x) {

                        register __m128 xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_ps(x.re,x.re);
                        xmm1 = _mm_mul_ps(x.im,x.im);
                        cabs = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        return (cabs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_xmm4c4_u(const float * __restrict re,
                                       const float * __restrict im,
                                       float * __restrict  carg) {

                        register __m128 xmm0,xmm1;
                        xmm0 = _mm_loadu_ps(&re[0]);
                        xmm1 = _mm_loadu_ps(&im[0]);
                        _mm_storeu_ps(&carg[0],_mm_atan2_ps(xmm0,xmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) re,
                                       const float * __restrict __ATTR_ALIGN__(16) im,
                                       float * __restrict  __ATTR_ALIGN__(16) carg) {

                        register __m128 xmm0,xmm1;
                        xmm0 = _mm_load_ps(&re[0]);
                        xmm1 = _mm_load_ps(&im[0]);
                        _mm_store_ps(&carg[0],_mm_atan2_ps(xmm0,xmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m128 carg_xmm4c4(const __m128 re,
                                       const __m128 im) {

                       register __m128 carg;
                       carg = _mm_atan2_ps(re,im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m128 carg_xmm4c4(xmm4c4_t x) {

                       register __m128 carg;
                       carg = _mm_atan2_ps(x.re,x.im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_xmm4c4(const __m128 re,
	                             const __m128 im,
	                             __m128 * __restrict clogr,
	                             __m128 * __restrict clogi) {
	                
	                register __m128 t1,t2,ln;
	                t1  = cabs_xmm4c4(re,im);
	                t2  = carg_xmm4c4(re,im);
	                ln  = _mm_log_ps(t1);
	                *clogr = ln;
	                *clogi = t2;                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) pre,
	                               const float * __restrict __ATTR_ALIGN__(16) pim,
	                               float * __restrict clogr,
	                               float * __restrict clogi) {
	                
	                register __m128 re = _mm_load_ps(&pre[0]);
	                register __m128 im = _mm_load_ps(&pim[0]);
	                register __m128 t1,t2,ln;
	                t1  = cabs_xmm4c4(re,im);
	                t2  = carg_xmm4c4(re,im);
	                ln  = _mm_log_ps(t1);
	                _mm_store_ps(&clogr[0] ,ln);
	                _mm_store_ps(&clogi[0] ,t2);                    
	        }
	        
	        
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_xmm4c4_u(const float * __restrict  pre,
	                               const float * __restrict  pim,
	                               float * __restrict clogr,
	                               float * __restrict clogi) {
	                
	                register __m128 re = _mm_loadu_ps(&pre[0]);
	                register __m128 im = _mm_loadu_ps(&pim[0]);
	                register __m128 t1,t2,ln;
	                t1  = cabs_xmm4c4(re,im);
	                t2  = carg_xmm4c4(re,im);
	                ln  = _mm_log_ps(t1);
	                _mm_storeu_ps(&clogr[0] ,ln);
	                _mm_storeu_ps(&clogi[0] ,t2);                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           xmm4c4_t clog_xmm4c4(const xmm4c4_t x){
	                                  
	                xmm4c4_t clog;                           
	                register __m128 t1,t2,ln;
	                t1  = cabs_xmm4c4(x.re,x.im);
	                t2  = carg_xmm4c4(x.re,x.im);
	                ln  = _mm_log_ps(t1);
	                clog.re = ln;
	                clog.im = t2;
	                return (clog);                    
	        }

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_xmm4c4_u(float * __restrict re,
                                        float * __restrict im) {

                        register __m128 c;
                        c = negate_xmm4r4(_mm_loadu_ps(&im[0]));
                        _mm_storeu_ps(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_xmm4c4_a(float * __restrict __ATTR_ALIGN__(16) re,
                                        float * __restrict __ATTR_ALIGN__(16) im) {
                                        
                        register __m128 c;
                        c = negate_xmm4r4(_mm_load_ps(&im[0]));
                        _mm_store_ps(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_xmm4c4(__m128 * __restrict re,
                                      __m128 * __restrict im) {
                         
                        register __m128 c;              
                        c = negate_xmm4r4(*im);
                        *im = c;
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_xmm4c4_v2(const __m128 xre,
                                         const __m128 xim,
                                         __m128 * __restrict yre,
                                         __m128 * __restrict yim) {
                         
                        //register __m128 c;              
                        //c = negate_xmm4c4(*im);
                        //*im = c;
                        *yre = xre; 
                        *yim = negate_xmm4r4(xim);
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cconj_xmm4c4_v2(const __m128 xre,
                                              const __m128 xim) {                                              
                         
                        //register __m128 c;              
                        //c = negate_xmm4c4(*im);
                        //*im = c;
                        xmm4c4_t cv;
                        cv.re = xre; 
                        cv.im = negate_xmm4r4(xim);
                        return (cv);
                   } 
                   
                   
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cconj_xmm4c4_v2(const xmm4c4_t x) {                                              
                         
                        //register __m128 c;              
                        //c = negate_xmm4c4(*im);
                        //*im = c;
                        xmm4c4_t cv;
                        cv.re = x.re; 
                        cv.im = negate_xmm4r4(x.im);
                        return (cv);
                   } 
                   
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict  csre,
                                       float * __restrict  csim) {

                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_loadu_ps(&xre[0]);
                      xmm1  = _mm_loadu_ps(&xim[0]);
                      xmm2  = _mm_mul_ps(_mm_cos_ps(xmm0),_mm_cosh_ps(xmm1));
                      _mm_storeu_ps(&csre[0],xmm2);
                      xmm3  = _mm_mul_ps(_mm_sin_ps(xmm0),_mm_sinh_ps(xmm1));
                      _mm_storeu_ps(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       float * __restrict  __ATTR_ALIGN__(16) csre,
                                       float * __restrict  __ATTR_ALIGN__(16) csim) {

                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_load_ps(&xre[0]);
                      xmm1  = _mm_load_ps(&xim[0]);
                      xmm2  = _mm_mul_ps(_mm_cos_ps(xmm0),_mm_cosh_ps(xmm1));
                      _mm_store_ps(&csre[0],xmm2);
                      xmm3  = _mm_mul_ps(_mm_sin_ps(xmm0),_mm_sinh_ps(xmm1));
                      _mm_store_ps(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     __m128 * __restrict csre,
                                     __m128 * __restrict csim) {

                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cos_ps(xre),_mm_cosh_ps(xim));
                      *csre = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sin_ps(xre),_mm_sinh_ps(xim));
                      *csim = xmm1; 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t ccos_xmm4c4(const __m128 xre,
                                          const __m128 xim) {
                                    
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cos_ps(xre),_mm_cosh_ps(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sin_ps(xre),_mm_sinh_ps(xim));
                      cv.im = xmm1;
                      return (cv); 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t ccos_xmm4c4(const xmm4c4_t x) {
                                    
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cos_ps(x.re),_mm_cosh_ps(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sin_ps(x.re),_mm_sinh_ps(x.im));
                      cv.im = xmm1;
                      return (cv); 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict  csre,
                                       float * __restrict  csim) {

                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_loadu_ps(&xre[0]);
                      xmm1  = _mm_loadu_ps(&xim[0]);
                      xmm2  = _mm_mul_ps(_mm_cosh_ps(xmm0),_mm_cos_ps(xmm1));
                      _mm_storeu_ps(&csre[0],xmm2);
                      xmm3  = _mm_mul_ps(_mm_sinh_ps(xmm0),_mm_sin_ps(xmm1));
                      _mm_storeu_ps(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       float * __restrict  __ATTR_ALIGN__(16) csre,
                                       float * __restrict  __ATTR_ALIGN__(16) csim) {

                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      xmm0  = _mm_load_ps(&xre[0]);
                      xmm1  = _mm_load_ps(&xim[0]);
                      xmm2  = _mm_mul_ps(_mm_cosh_ps(xmm0),_mm_cos_ps(xmm1));
                      _mm_store_ps(&csre[0],xmm2);
                      xmm3  = _mm_mul_ps(_mm_sinh_ps(xmm0),_mm_sin_ps(xmm1));
                      _mm_store_ps(&csim[0],xmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     __m128 * __restrict csre,
                                     __m128 * __restrict csim) {

                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cosh_ps(xre),_mm_cos_ps(xim));
                      *csre = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sinh_ps(xre),_mm_sin_ps(xim));
                      *csim = xmm1; 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t ccosh_xmm4c4(const __m128 xre,
                                           const __m128 xim) {
                                          
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cosh_ps(xre),_mm_cos_ps(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sinh_ps(xre),_mm_sin_ps(xim));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t ccosh_xmm4c4(const xmm4c4_t x) {
                                          
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cosh_ps(x.re),_mm_cos_ps(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sinh_ps(x.re),_mm_sin_ps(x.im));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_xmm4c4(const __m128 xre,
	                             const __m128 xim,
	                             const float n,
	                             __m128 * __restrict powr,
	                             __m128 * __restrict powi) {
	                             
	                register __m128 xmm0,xmm1;
	                register __m128 r,tht;
	                register __m128 vn,pt;
	                register __m128 ta;
	                xmm0  = _mm_mul_ps(xre,xre);
	                vn    = _mm_set1_ps(n);
	                xmm1  = _mm_mul_ps(xim,xim);
	                r     = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
	                tht   = _mm_atan_ps(_mm_div_ps(xim,xre));
	                pt    = _mm_pow_ps(r,vn);
	                ta    = _mm_mul_ps(vn,tht);
	                *powr = _mm_mul_ps(pt,_mm_cos_ps(ta));
	                *powi = _mm_mul_ps(pt,_mm_sin_ps(ta));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) pxre,
	                               const float * __restrict __ATTR_ALIGN__(16) pxim,
	                               const float n,
	                               float * __restrict __ATTR_ALIGN__(16) ppowr
	                               float * __restrict __ATTR_ALIGN__(16) ppowi) {
	                  
	                register __m128 xre = _mm_load_ps(&pxre[0]);
	                register __m128 xim = _mm_load_ps(&pxim[0]);           
	                register __m128 xmm0,xmm1;
	                register __m128 r,tht;
	                register __m128 vn,pt;
	                register __m128 ta;
	                xmm0  = _mm_mul_ps(xre,xre);
	                vn    = _mm_set1_ps(n);
	                xmm1  = _mm_mul_ps(xim,xim);
	                r     = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
	                tht   = _mm_atan_ps(_mm_div_ps(xim,xre));
	                pt    = _mm_pow_ps(r,vn);
	                ta    = _mm_mul_ps(vn,tht);
	                _mm_store_ps(&ppowr[0] ,_mm_mul_ps(pt,_mm_cos_ps(ta)));
	                _mm_store_ps(&ppowi[0] ,_mm_mul_ps(pt,_mm_sin_ps(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_xmm4c4_u(const float * __restrict  pxre,
	                               const float * __restrict  pxim,
	                               const float n,
	                               float * __restrict  ppowr
	                               float * __restrict  ppowi) {
	                  
	                register __m128 xre = _mm_loadu_ps(&pxre[0]);
	                register __m128 xim = _mm_loadu_ps(&pxim[0]);           
	                register __m128 xmm0,xmm1;
	                register __m128 r,tht;
	                register __m128 vn,pt;
	                register __m128 ta;
	                xmm0  = _mm_mul_ps(xre,xre);
	                vn    = _mm_set1_ps(n);
	                xmm1  = _mm_mul_ps(xim,xim);
	                r     = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
	                tht   = _mm_atan_ps(_mm_div_ps(xim,xre));
	                pt    = _mm_pow_ps(r,vn);
	                ta    = _mm_mul_ps(vn,tht);
	                _mm_storeu_ps(&ppowr[0] ,_mm_mul_ps(pt,_mm_cos_ps(ta)));
	                _mm_storeu_ps(&ppowi[0] ,_mm_mul_ps(pt,_mm_sin_ps(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ZMM16c4_t cpow_xmm4c4(const ZMM16c4_t x,
	                                  const float n) {
	                   
	                ZMM16c4_t cp;        
	                register __m128 xmm0,xmm1;
	                register __m128 r,tht;
	                register __m128 vn,pt;
	                register __m128 ta;
	                xmm0  = _mm_mul_ps(x.re,x.re);
	                vn    = _mm_set1_ps(n);
	                xmm1  = _mm_mul_ps(x.im,x.im);
	                r     = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
	                tht   = _mm_atan_ps(_mm_div_ps(x.im,x.re));
	                pt    = _mm_pow_ps(r,vn);
	                ta    = _mm_mul_ps(vn,tht);
	                cp.re = _mm_mul_ps(pt,_mm_cos_ps(ta));
	                cp.im = _mm_mul_ps(pt,_mm_sin_ps(ta));      
	                return (cp);              
	       }
	       
#include <utility>

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm4c4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim) {
                                      
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_loadu_ps(&xre[0]);
                      xmm1 = _mm_loadu_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_EQ_OQ);
                      xmm2 = _mm_loadu_ps(&xim[0]);
                      xmm3 = _mm_loadu_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_EQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                      const float * __restrict __ATTR_ALIGN__(16) xim,
                                      const float * __restrict __ATTR_ALIGN__(16) yre,
                                      const float * __restrict __ATTR_ALIGN__(16) yim) {
                                      
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_load_ps(&xre[0]);
                      xmm1 = _mm_load_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_EQ_OQ);
                      xmm2 = _mm_load_ps(&xim[0]);
                      xmm3 = _mm_load_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_EQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm4c4(      const __m128 xre,
                                    const __m128 xim,
                                    const __m128 yre,
                                    const __m128 yim) {
                                    
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm_cmp_ps_mask(xre,yre,_CMP_EQ_OQ);
                         eqi = _mm_cmp_ps_mask(xim,yim,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm4c4(      const xmm4c4_t x,
                                    const xmm4c4_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm4c4_u(      const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim) {
                                     
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_loadu_ps(&xre[0]);
                      xmm1 = _mm_loadu_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_GT_OQ);
                      xmm2 = _mm_loadu_ps(&xim[0]);
                      xmm3 = _mm_loadu_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_GT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm4c4_a(      const float * __restrict __ATTR_ALIGN__(16) xre,
                                      const float * __restrict __ATTR_ALIGN__(16) xim,
                                      const float * __restrict __ATTR_ALIGN__(16) yre,
                                      const float * __restrict __ATTR_ALIGN__(16) yim) {
                                      
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_load_ps(&xre[0]);
                      xmm1 = _mm_load_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_GT_OQ);
                      xmm2 = _mm_load_ps(&xim[0]);
                      xmm3 = _mm_load_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_GT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm4c4(      const __m128 xre,
                                    const __m128 xim,
                                    const __m128 yre,
                                    const __m128 yim) {
                            
                         __mmask8 eqr;
                         __mmask8 eqi;      
                         eqr = _mm_cmp_ps_mask(xre,yre,_CMP_GT_OQ);
                         eqi = _mm_cmp_ps_mask(xim,yim,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm4c4(      const xmm4c4_t x,
                                    const xmm4c4_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm4c4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim){
                                      
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      xmm0 = _mm_loadu_ps(&xre[0]);
                      xmm1 = _mm_loadu_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_LT_OQ);
                      xmm2 = _mm_loadu_ps(&xim[0]);
                      xmm3 = _mm_loadu_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_LT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm4c4_a(      const float * __restrict __ATTR_ALIGN__(16) xre,
                                      const float * __restrict __ATTR_ALIGN__(16) xim,
                                      const float * __restrict __ATTR_ALIGN__(16) yre,
                                      const float * __restrict __ATTR_ALIGN__(16) yim) {
                                     
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      xmm0 = _mm_load_ps(&xre[0]);
                      xmm1 = _mm_load_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_LT_OQ);
                      xmm2 = _mm_load_ps(&xim[0]);
                      xmm3 = _mm_load_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_LT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm4c4(      const __m128 xre,
                                    const __m128 xim,
                                    const __m128 yre,
                                    const __m128 yim) {
                               
                         __mmask8 eqr;
                         __mmask8 eqi;      
                         eqr = _mm_cmp_ps_mask(xre,yre,_CMP_LT_OQ);
                         eqi = _mm_cmp_ps_mask(xim,yim,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm4c4(     const xmm4c4_t x,
                                    const xmm4c4_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm4c4_u(     const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim) {
                                      
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      xmm0 = _mm_loadu_ps(&xre[0]);
                      xmm1 = _mm_loadu_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_NEQ_OQ);
                      xmm2 = _mm_loadu_ps(&xim[0]);
                      xmm3 = _mm_loadu_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_NEQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                      const float * __restrict __ATTR_ALIGN__(16) xim,
                                      const float * __restrict __ATTR_ALIGN__(16) yre,
                                      const float * __restrict __ATTR_ALIGN__(16) yim) {
                                      
                      register __m128 xmm0,xmm1,xmm2,xmm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      xmm0 = _mm_load_ps(&xre[0]);
                      xmm1 = _mm_load_ps(&yre[0]);
                      eqr  = _mm_cmp_ps_mask(xmm0,xmm1,_CMP_NEQ_OQ);
                      xmm2 = _mm_load_ps(&xim[0]);
                      xmm3 = _mm_load_ps(&yim[0]);
                      eqi  = _mm_cmp_ps_mask(xmm2,xmm3,_CMP_NEQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm4c4(     const __m128 xre,
                                    const __m128 xim,
                                    const __m128 yre,
                                    const __m128 yim) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm_cmp_ps_mask(xre,yre,_CMP_NEQ_OQ);
                         eqi = _mm_cmp_ps_mask(xim,yim,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm4c4(      const xmm4c4_t x,
                                     const xmm4c4_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict cexpr,
                                       float * __restrict cexpi ) {

                        register const __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_loadu_ps(&xre[0]);
                        xmm1  = _mm_loadu_ps(&xim[0]);
                        xmm2  = _mm_exp_ps(xmm0);
                        xmm3  = _mm_mul_ps(xmm2,_mm_cos_ps(xmm1));
                        _mm_storeu_ps(&cexpr[0],xmm3);
                        xmm4  = _mm_mul_ps(xmm2,_mm_mul_ps(_mm_sin_ps(xmm1),I));
                        _mm_storeu_ps(&cexpi[0],xmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       float * __restrict __ATTR_ALIGN__(16) cexpr,
                                       float * __restrict __ATTR_ALIGN__(16) cexpi ) {

                        register const __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        xmm0  = _mm_load_ps(&xre[0]);
                        xmm1  = _mm_load_ps(&xim[0]);
                        xmm2  = _mm_exp_ps(xmm0);
                        xmm3  = _mm_mul_ps(xmm2,_mm_cos_ps(xmm1));
                        _mm_store_ps(&cexpr[0],xmm3);
                        xmm4  = _mm_mul_ps(xmm2,_mm_mul_ps(_mm_sin_ps(xmm1),I));
                        _mm_store_ps(&cexpi[0],xmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_xmm4c4(const __m128 xre,
                                     const __m128 xim,
                                     __m128 * __restrict cexpr,
                                     __m128 * __restrict cexpi) {

                        register const __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0;
                        xmm0   = _mm_exp_ps(xre);
                        *cexpr = _mm_mul_ps(xmm0,_mm_cos_ps(xim));
                        *cexpi = _mm_mul_ps(xmm0,_mm_mul_ps(_mm_sin_ps(xim),I));
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cexp_xmm4c4(const __m128 xre,
                                          const __m128 xim) {
                                     
                        xmm4c4_t cv;
                        register const __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0;
                        xmm0   = _mm_exp_ps(xre);
                        cv.re = _mm_mul_ps(xmm0,_mm_cos_ps(xim));
                        cv.im = _mm_mul_ps(xmm0,_mm_mul_ps(_mm_sin_ps(xim),I));
                        return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cexp_xmm4c4(const xmm4c4_t x) {
                                     
                        xmm4c4_t cv;
                        register const __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0;
                        xmm0   = _mm_exp_ps(x.re);
                        cv.re = _mm_mul_ps(xmm0,_mm_cos_ps(x.im));
                        cv.im = _mm_mul_ps(xmm0,_mm_mul_ps(_mm_sin_ps(x.im),I));
                        return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_xmm4c4_u(const float * __restrict rho,
                                         const float * __restrict tht,
                                         float * __restrict  re,
                                         float * __restrict  im) {

                         register __m128 xmm0,xmm1,xmm2,xmm3;
                         xmm0 = _mm_loadu_ps(&rho[0]);
                         xmm1 = _mm_loadu_ps(&tht[0]);
                         xmm2 = _mm_mul_ps(xmm0,_mm_cos_ps(xmm1)); //tht
                         _mm_storeu_ps(&re[0],xmm2);
                         xmm3 = _mm_mul_ps(xmm0,_mm_sin_ps(xmm1)); //tht
                         _mm_storeu_ps(&im[0],xmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) rho,
                                         const float * __restrict __ATTR_ALIGN__(16) tht,
                                         float * __restrict  __ATTR_ALIGN__(16) re,
                                         float * __restrict  __ATTR_ALIGN__(16) im) {

                         register __m128 xmm0,xmm1,xmm2,xmm3;
                         xmm0 = _mm_load_ps(&rho[0]);
                         xmm1 = _mm_load_ps(&tht[0]);
                         xmm2 = _mm_mul_ps(xmm0,_mm_cos_ps(xmm1)); //tht
                         _mm_store_ps(&re[0],xmm2);
                         xmm3 = _mm_mul_ps(xmm0,_mm_sin_ps(xmm1)); //tht
                         _mm_store_ps(&im[0],xmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_xmm4c4(const __m128 rho,
                                       const __m128 tht,
                                       __m128 * __restrict re,
                                       __m128 * __restrict im) {

                        register __m128 xmm0,xmm1;
                        xmm0 = _mm_mul_ps(rho,_mm_cos_ps(tht));
                        *re  = xmm0;
                        xmm1 = _mm_mul_ps(rho,_mm_sin_ps(tht));
                        *im  = xmm1;
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cpolar_xmm4c4(const __m128 rho,
                                            const __m128 tht) {
                                      
                        xmm4c4_t cv
                        register __m128 xmm0,xmm1;
                        xmm0 = _mm_mul_ps(rho,_mm_cos_ps(tht));
                        cv.re  = xmm0;
                        xmm1 = _mm_mul_ps(rho,_mm_sin_ps(tht));
                        cv.im  = xmm1;
                        return (cv);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_xmm4c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict wrkc,
                                       float * __restrict csqr,
                                       float * __restrict csqi) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 half = _mm_set1_ps(0.5);
                        cabs_xmm4c4_u(xre,xim,wrkc);
                        xmm0  = _mm_loadu_ps(&xre[0]);
                        xmm1  = _mm_loadu_ps(&wrkc[0]);
                        xmm2  = _mm_mul_ps(half,_mm_add_ps(xmm1,xmm0));
                        _mm_storeu_ps(&csqr[0],_mm_sqrt_ps(xmm2));
                        xmm3  = _mm_mul_ps(half,_mm_sub_ps(xmm1,xmm0));
                        _mm_storeu_ps(&csqi[0],_mm_sqrt_ps(xmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                       const float * __restrict __ATTR_ALIGN__(16) xim,
                                       float * __restrict __ATTR_ALIGN__(16) wrkc,
                                       float * __restrict __ATTR_ALIGN__(16) csqr,
                                       float * __restrict __ATTR_ALIGN__(16) csqi) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 half = _mm_set1_ps(0.5);
                        cabs_xmm4c4_a(xre,xim,wrkc);
                        xmm0  = _mm_load_ps(&xre[0]);
                        xmm1  = _mm_load_ps(&wrkc[0]);
                        xmm2  = _mm_mul_ps(half,_mm_add_ps(xmm1,xmm0));
                        _mm_store_ps(&csqr[0],_mm_sqrt_ps(xmm2));
                        xmm3  = _mm_mul_ps(half,_mm_sub_ps(xmm1,xmm0));
                        _mm_store_ps(&csqi[0],_mm_sqrt_ps(xmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_xmm4c4(const __m128 xre,
                                      const __m128 xim,
                                      __m128 * __restrict wrkc,
                                      __m128 * __restrict csqr,
                                      __m128 * __restrict csqi) {

                       register __m128 xmm0,xmm1;
                       register __m128 half = _mm_set1_ps(0.5); 
                       cabs_xmm4c4(xre,xim,wrkc);
                       xmm0  = _mm_mul_ps(half,_mm_add_ps(*wrkc,xre));
                       *csqr = xmm0;
                       xmm1  = _mm_mul_ps(half,_mm_sub_ps(*wrkc,xre));
                       *csqi = xmm1; 
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t csqrt_xmm4c4(const __m128 xre,
                                           const __m128 xim,
                                          __m128 * __restrict wrkc) {
                                          
                       xmm4c4_t cv;
                       register __m128 xmm0,xmm1;
                       register __m128 half = _mm_set1_ps(0.5); 
                       cabs_xmm4c4(xre,xim,wrkc);
                       xmm0  = _mm_mul_ps(half,_mm_add_ps(*wrkc,xre));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_ps(half,_mm_sub_ps(*wrkc,xre));
                       cv.im = xmm1; 
                       return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t csqrt_xmm4c4(const xmm4c4_t x,
                                          __m128 * __restrict wrkc) {
                                          
                       xmm4c4_t cv;
                       register __m128 xmm0,xmm1;
                       register __m128 half = _mm_set1_ps(0.5); 
                       cabs_xmm4c4(x.re,x.im,wrkc);
                       xmm0  = _mm_mul_ps(half,_mm_add_ps(*wrkc,x.re));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_ps(half,_mm_sub_ps(*wrkc,x.re));
                       cv.im = xmm1; 
                       return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_xmm4c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict zre,
                                             float * __restrict zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        rep  = _mm_fmsub_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        imp  = _mm_fmadd_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        _mm_storeu_ps(&zre[0], _mm_div_ps(rep,xmm2));
                        _mm_storeu_ps(&zim[0], _mm_div_ps(imp,xmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                             const float * __restrict __ATTR_ALIGN__(16) xim,
                                             const float * __restrict __ATTR_ALIGN__(16) yre,
                                             const float * __restrict __ATTR_ALIGN__(16) yim,
                                             float * __restrict __ATTR_ALIGN__(16) zre,
                                             float * __restrict __ATTR_ALIGN__(16) zim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        rep  = _mm_fmsub_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        imp  = _mm_fmadd_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        _mm_store_ps(&zre[0], _mm_div_ps(rep,xmm2));
                        _mm_store_ps(&zim[0], _mm_div_ps(imp,xmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_xmm4c4(  const __m128  xre,
                                             const __m128  xim,
                                             const __m128  yre,
                                             const __m128  yim,
                                             __m128 * __restrict zre,
                                             __m128 * __restrict zim) {

                        register __m128 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        imp  = _mm_fmadd_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        *zre = _mm_div_ps(rep,xmm2);
                        *zim = _mm_div_ps(imp,xmm2);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cnorm_prod_xmm4c4(  const __m128  xre,
                                                  const __m128  xim,
                                                  const __m128  yre,
                                                  const __m128  yim) {
                                             
                        xmm4c4_t cv;
                        register __m128 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        imp  = _mm_fmadd_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        cv.re = _mm_div_ps(rep,xmm2);
                        cv.im = _mm_div_ps(imp,xmm2);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cnorm_prod_xmm4c4(  const xmm4c4_t x,
                                                  const xmm4c4_t y) {
                                             
                        xmm4c4_t cv;
                        register __m128 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        imp  = _mm_fmadd_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        cv.re = _mm_div_ps(rep,xmm2);
                        cv.im = _mm_div_ps(imp,xmm2);
                        return (cv);
             }
             



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_xmm4c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        sre = 0.0;
                        rep  = _mm_fmsub_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        t0   = _mm_hadd_ps(rep,rep);//xmm8r4_horizontal_sum(rep);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        t1   = _mm_hadd_ps(imp,imp);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_xmm4c4_a( const float * __restrict __ATTR_ALIGN__(16) xre,
                                             const float * __restrict __ATTR_ALIGN__(16) xim,
                                             const float * __restrict __ATTR_ALIGN__(16) yre,
                                             const float * __restrict __ATTR_ALIGN__(16) yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        sre = 0.0;
                        rep  = _mm_fmsub_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        t0   = _mm_hadd_ps(rep,rep);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        t1   = _mm_hadd_ps(imp,imp);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
              } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_xmm4c4( const __m128 xre,
                                           const __m128 xim,
                                           const __m128 yre,
                                           const __m128 yim,
                                           float * __restrict mre,
                                           float * __restrict min) {

                        register __m128 rep,imp;
                        __m128 t0,t1;
                        constexpr float inv16 = 0.5;
                        float sre,sim;
                        sre = 0.0;
                        rep  = _mm_fmsub_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        t0   = _mm_hadd_ps(rep,rep);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        t1   = _mm_hadd_ps(imp,imp);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_xmm4c4(const xmm4c4_t x,
                                           const xmm4c4_t y,
                                           float * __restrict mre,
                                           float * __restrict min) {

                        register __m128 rep,imp;
                        __m128 t0,t1
                        constexpr float inv16 = 0.5;
                        float sre,sim;
                        sre = 0.0;
                        rep  = _mm_fmsub_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        t0   = _mm_hadd_ps(rep,rep);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        sim  = 0.0;
                        imp  = _mm_fmadd_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        t1   = _mm_hadd_ps(imp,imp);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }


               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_xmm4c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp,den,rquot,iquot;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        sre  = 0.0;
                        rep  = _mm_fmsub_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        imp  = _mm_fmadd_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        sim  = 0.0;
                        den  = _mm_fmadd_ps(xmm1,xmm1,
                                               _mm_mul_ps(xmm3,xmm3));
                        rquot = _mm_div_ps(rep,den);
                        t0    = _mm_hadd_ps(rquot,rquot);
                        sre   = *(float*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_ps(imp,den);
                        t1    = _mm_hadd_ps(iquot,iquot);
                        sim   = *(float*)&t1;
                        *mim  = sre*inv2;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                             const float * __restrict __ATTR_ALIGN__(16) xim,
                                             const float * __restrict __ATTR_ALIGN__(16) yre,
                                             const float * __restrict __ATTR_ALIGN__(16) yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp,den,rquot,iquot;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm_fmsub_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        imp  = _mm_fmadd_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        sim  = 0.0f;
                        den  = _mm_fmadd_ps(xmm1,xmm1,
                                               _mm_mul_ps(xmm3,xmm3));
                        rquot = _mm_div_ps(rep,den);
                        t0    = _mm_hadd_ps(rquot,rquot);
                        sre   = *(float*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_ps(imp,den);
                        t1    = _mm_hadd_ps(iquot,iquot);
                        sim   = *(float*)&t1;
                        *mim  = sre*inv2;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_xmm4c4(  const __m128 xre,
                                             const __m128 xim,
                                             const __m128 yre,
                                             const __m128 yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 rep,imp,den,rquot,iquot;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        sre  = 0.0;
                        rep  = _mm_fmsub_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        imp  = _mm_fmadd_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        sim  = 0.0;
                        den  = _mm_fmadd_ps(yre,yre,
                                               _mm_mul_ps(yim,yim));
                        rquot = _mm_div_ps(rep,den);
                        t0    = _mm_hadd_ps(rquot,rquot);
                        sre   = *(float*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_ps(imp,den);
                        t1    = _mm_hadd_ps(iquot,iquot);
                        sim   = *(float*)&t1;
                        *mim  = sre*inv2;
              }  
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_xmm4c4(  const xmm4c4_t x,
                                             const xmm4c4_t y,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 rep,imp,den,rquot,iquot;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        sre  = 0.0;
                        rep  = _mm_fmsub_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        imp  = _mm_fmadd_ps(xim,yre,
                                               _mm_mul_ps(x.re,y.im));
                        sim  = 0.0;
                        den  = _mm_fmadd_ps(y.re,y.re,
                                               _mm_mul_ps(y.im,y.im));
                        rquot = _mm_div_ps(rep,den);
                        t0    = _mm_hadd_ps(rquot,rquot);
                        sre   = *(float*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_ps(imp,den);
                        t1    = _mm_hadd_ps(iquot,iquot);
                        sim   = *(float*)&t1;
                        *mim  = sre*inv2;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_xmm4c4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim ) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp,magc1,magc2,vcmag;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        rep  = _mm_fmadd_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        _mm_storeu_ps(&mre[0], _mm_div_ps(rep,vcmag));
                        _mm_storeu_ps(&mim[0], _mm_div_ps(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                              const float * __restrict __ATTR_ALIGN__(16) xim,
                                              const float * __restrict __ATTR_ALIGN__(16) yre,
                                              const float * __restrict __ATTR_ALIGN__(16) yim,
                                              float * __restrict __ATTR_ALIGN__(16) mre,
                                              float * __restrict __ATTR_ALIGN__(16) mim ) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 rep,imp,magc1,magc2,vcmag;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        rep  = _mm_fmad_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        _mm_store_ps(&mre[0], _mm_div_ps(rep,vcmag));
                        _mm_store_ps(&mim[0], _mm_div_ps(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_xmm4c4(const __m128 xre,
                                            const __m128 xim,
                                            const __m128 yre,
                                            const __m128 yim,
                                            __m128 * __restrict mre,
                                            __m128 * __restrict mim) {

                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        *mre = _mm_div_ps(rep,vcmag);
                        *mim = _mm_div_ps(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_xmm4c4(const xmm4c4_t x,
                                            const xmm4c4_t y,
                                            __m128 * __restrict mre,
                                            __m128 * __restrict mim) {

                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        *mre = _mm_div_ps(rep,vcmag);
                        *mim = _mm_div_ps(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cnorm_cprod_xmm4c4(const __m128 xre,
                                                 const __m128 xim,
                                                 const __m128 yre,
                                                 const __m128 yim) {
                                               
                        xmm4c4_t cv;
                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        cv.re = _mm_div_ps(rep,vcmag);
                        cv.im = _mm_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cnorm_cprod_xmm4c4(const xmm4c4_t x,
                                                 const xmm4c4_t y) {
                                               
                        xmm4c4_t cv;
                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmad_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        cv.re = _mm_div_ps(rep,vcmag);
                        cv.im = _mm_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
             


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_xmm4c4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim) {
                      
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        re   = _mm_fmadd_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                              const float * __restrict __ATTR_ALIGN__(16) xim,
                                              const float * __restrict __ATTR_ALIGN__(16) yre,
                                              const float * __restrict __ATTR_ALIGN__(16) yim,
                                              float * __restrict mre,
                                              float * __restrict mim) {
                      
                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        re   = _mm_fmadd_ps(xmm0,xmm1,
                                               _mm_mul_ps(xmm2,xmm3));
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_ps(xmm2,xmm1,
                                               _mm_mul_ps(xmm0,xmm3));
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_xmm4c4( const __m128 xre,
                                            const __m128 xim,
                                            const __m128 yre,
                                            const __m128 yim,
                                            float * __restrict mre,
                                            float * __restrict min) {

                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        re   = _mm_fmadd_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_xmm4c4(const xmm4c4_t x,
                                            const xmm4c4_t y,
                                            float * __restrict mre,
                                            float * __restrict min) {

                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.5;
                        float sre,sim;
                        re   = _mm_fmadd_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_xmm4c4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv16 = 0.5;
                        float sre,sim;
                        re   = _mm_loadu_ps(&xre[0]);
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_loadu_ps(&xim[0]);
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) xre,
                                              const float * __restrict __ATTR_ALIGN__(16) xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv16 = 0.5;
                        float sre,sim;
                        re   = _mm_load_ps(&xre[0]);
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_load_ps(&xim[0]);
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_xmm4c4(  const __m128 xre,
                                              const __m128 xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        constexpr float inv2 = 0.5;
                        __m128 t0,t1;
                        float sre,sim;
                        t0   = _mm_hadd_ps(xre,xre);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_ps(xim,xim);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_xmm4c4(  const xmm4c4_t x,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        constexpr float inv2 = 0.5;
                        __m128 t0,t1;
                        float sre,sim;
                        t0   = _mm_hadd_ps(x.re,x.re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_ps(x.im,x.im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_xmm4c4_u( const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim ) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 re,im,cvmag;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xmm0,xmm1,
                                                              _mm_mul_ps(xmm2,xmm3)));
                        _mm_storeu_ps(&mre[0], _mm_div_ps(xmm0,cvmag));
                        _mm_storeu_ps(&mim[0], _mm_div_ps(xmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_xmm4c4_a( const float * __restrict __ATTR_ALIGN__(16) xre,
                                              const float * __restrict __ATTR_ALIGN__(16) xim,
                                              const float * __restrict __ATTR_ALIGN__(16) yre,
                                              const float * __restrict __ATTR_ALIGN__(16) yim,
                                              float * __restrict __ATTR_ALIGN__(16) mre,
                                              float * __restrict __ATTR_ALIGN__(16) mim ) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 re,im,cvmag;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xmm0,xmm1,
                                                              _mm_mul_ps(xmm2,xmm3)));
                        _mm_store_ps(&mre[0], _mm_div_ps(xmm0,cvmag));
                        _mm_store_ps(&mim[0], _mm_div_ps(xmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_xmm4c4( const __m128 xre,
                                            const __m128 xim,
                                            const __m128 yre,
                                            const __m128 yim,
                                            __m128 * __restrict mre,
                                            __m128 * __restrict mim ) {

                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xre,yre,
                                                    _mm_mul_ps(xim,yim)));
                        *mre = _mm_div_ps(xre,cvmag));
                        *mim =  _mm_div_ps(xim,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_xmm4c4( const xmm4c4_t x,
                                            const xmm4c4_t y,
                                            __m128 * __restrict mre,
                                            __m128 * __restrict mim ) {

                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(x.re,y.re,
                                                    _mm_mul_ps(x.im,y.im)));
                        *mre = _mm_div_ps(x.re,cvmag));
                        *mim =  _mm_div_ps(x.im,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cnormalize_xmm4c4( const __m128 xre,
                                                 const __m128 xim,
                                                 const __m128 yre,
                                                 const __m128 yim) {
                                            
                        xmm4c4_t cv;
                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xre,yre,
                                                    _mm_mul_ps(xim,yim)));
                        cv.re = _mm_div_ps(xre,cvmag));
                        cv.im =  _mm_div_ps(xim,cvmag));
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   xmm4c4_t cnormalize_xmm4c4( const xmm4c4_t x,
                                                 const xmm4c4_t y,) {
                                            
                        xmm4c4_t cv;
                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(x.re,y.re,
                                                    _mm_mul_ps(x.im,y.im)));
                        cv.re = _mm_div_ps(x.re,cvmag));
                        cv.im =  _mm_div_ps(x.im,cvmag));
                        return (cv);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_xmm4c4_u( const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 cvmag;
                        xmm0 = _mm_loadu_ps(&xre[0]);
                        xmm1 = _mm_loadu_ps(&yre[0]);
                        xmm2 = _mm_loadu_ps(&xim[0]);
                        xmm3 = _mm_loadu_ps(&yim[0]);
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xmm0,xmm1,
                                                          _mm_mul_ps(xmm2,xmm3)));
                        _mm_storeu_ps(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_xmm4c4_a( const float * __restrict __ATTR_ALIGN__(16) xre,
                                              const float * __restrict __ATTR_ALIGN__(16) xim,
                                              const float * __restrict __ATTR_ALIGN__(16) yre,
                                              const float * __restrict __ATTR_ALIGN__(16) yim,
                                              float * __restrict __ATTR_ALIGN__(16) mre) {

                        register __m128 xmm0,xmm1,xmm2,xmm3;
                        register __m128 cvmag;
                        xmm0 = _mm_load_ps(&xre[0]);
                        xmm1 = _mm_load_ps(&yre[0]);
                        xmm2 = _mm_load_ps(&xim[0]);
                        xmm3 = _mm_load_ps(&yim[0]);
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xmm0,xmm1,
                                                          _mm_mul_ps(xmm2,xmm3)));
                        _mm_store_ps(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_xmm4c4(   const __m128 xre,
                                              const __m128 xim,
                                              const __m128 yre,
                                              const __m128 yim,
                                              __m128 * __restrict  mre) {

                        register __m128 cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xre,yre,
                                                          _mm_mul_ps(xim,yim)));
                        *mre = cvmag;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_xmm4c4(   const xmm4c4_t x,
                                              const xmm4c4_t y,
                                              __m128 * __restrict  mre) {

                        register __m128 cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(x.re,y.re,
                                                          _mm_mul_ps(x.im,y.im)));
                        *mre = cvmag;
             }


#include <complex>


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_2xr4_c4_unroll16x(float * __restrict __ATTR_ALIGN__(64) xre,
                                               float * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<float> * __restrict __ATTR_ALIGN__(64) vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register float re0,im0,re1,im1,re2,im2;
                        register float re3,im3,re4,im4,re5,im5;
                        register float re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (float*)__builtin_assume_aligned(xre,64);
                        xim = (float*)__builtin_assume_aligned(xim,);
                        vc   = (std::complex<float>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const float re = xre[j];
                               const float im = xim[j];
                               const std::complex<float> tj{re,im};
                               vc[i] = tj;
                           }
                            
                            if(n<16) {return;}
                       }

                       m1 = m+1;
#if defined(__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(vc:64)
                       for(j = m1; j != n; j += 16) {
                           re0     = xre[i+0];
                           im0     = xim[i+0];
                           vc[i+0] = {re0,im0};
                           re1     = xre[i+1];
                           im1     = xim[i+1];
                           vc[i+1] = {re1,im1}; 
                           re2     = xre[i+2];
                           im2     = xim[i+2];
                           vc[i+2] = {re2,im2};
                           re3     = xre[i+3];
                           im3     = xim[i+3];
                           vc[i+3] = {re3,im3};
                           re4     = xre[i+4];
                           im4     = xim[i+4];
                           vc[i+4] = {re4,im4};
                           re5     = xre[i+5];
                           im5     = xim[i+5];
                           vc[i+5] = {re5,im5}; 
                           re6     = xre[i+6];
                           im6     = xim[i+6];
                           vc[i+6] = {re6,im6};
                           re7     = xre[i+7];
                           im7     = xim[i+7];
                           vc[i+7] = {re7,im7};
                           re0     = xre[i+8];
                           im0     = xim[i+8];
                           vc[i+8] = {re0,im0};
                           re1     = xre[i+9];
                           im1     = xim[i+9];
                           vc[i+9] = {re1,im1}; 
                           re2     = xre[i+10];
                           im2     = xim[i+10];
                           vc[i+10]= {re2,im2};
                           re3     = xre[i+11];
                           im3     = xim[i+11];
                           vc[i+11]= {re3,im3};
                           re4     = xre[i+12];
                           im4     = xim[i+12];
                           vc[i+12]= {re4,im4};
                           re5     = xre[i+13];
                           im5     = xim[i+13];
                           vc[i+13]= {re5,im5}; 
                           re6     = xre[i+14];
                           im6     = xim[i+14];
                           vc[i+14]= {re6,im6};
                           re7     = xre[i+15];
                           im7     = xim[i+15];
                           vc[i+15]= {re7,im7};
                      }
                 }


              

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_c4_2xr4_unroll16x( float * __restrict __ATTR_ALIGN__(64) xre,
                                               float * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<float> * __restrict __ATTR_ALIGN__(64) vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register float re0,im0,re1,im1,re2,im2;
                        register float re3,im3,re4,im4,re5,im5;
                        register float re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (float*)__builtin_assume_aligned(xre,64);
                        xim = (float*)__builtin_assume_aligned(xim,64);
                        vc   = (std::complex<float>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const float re = vc[j].real();
                               const float im = vc[j].imag();
                               xre[j]         = re;
                               xim[j]         = im;
                           }
                            
                            if(n<16) {return;}
                       }

                       m1 = m+1;
#if defined(__INTEL_COMPILER)
#pragma code_align(32)
#endif
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(vc:64)
                       for(j = m1; j != n; j += 16) {
                           re0    = vc[j+0].real();
                           xre[j+0] = re0;
                           im0    = vc[j+0].imag();
                           xim[j+0] = im0;
                           re1    = vc[j+1].real();
                           xre[j+1] = re1;
                           im1    = vc[j+1].imag();
                           xim[j+1] = im1;
                           re2    = vc[j+2].real();
                           xre[j+2] = re2;
                           im2    = vc[j+2].imag();
                           xim[j+2] = im2;
                           re3    = vc[j+3].real();
                           xre[j+3] = re3;
                           im3    = vc[j+3].imag();
                           xim[j+3] = im3;
                           re4    = vc[j+4].real();
                           xre[j+4] = re4;
                           im4    = vc[j+4].imag();
                           xim[j+4] = im4;
                           re5    = vc[j+5].real();
                           xre[j+5] = re5;
                           im5    = vc[j+5].imag();
                           xim[j+5] = im5;
                           re6    = vc[j+6].real();
                           xre[j+6] = re6;
                           im6    = vc[j+6].imag();
                           xim[j+6] = im6;
                           re7    = vc[j+7].real();
                           xre[j+7] = re7;
                           im7    = vc[j+7].imag();
                           xim[j+7] = im7;
                           re0    = vc[j+8].real();
                           xre[j+8] = re0;
                           im0    = vc[j+8].imag();
                           xim[j+8] = im0;
                           re1    = vc[j+9].real();
                           xre[j+9] = re1;
                           im1    = vc[j+9].imag();
                           xim[j+9] = im1;
                           re2    = vc[j+10].real();
                           xre[j+10] = re2;
                           im2    = vc[j+10].imag();
                           xim[j+10] = im2;
                           re3    = vc[j+11].real();
                           xre[j+11] = re3;
                           im3    = vc[j+11].imag();
                           xim[j+11] = im3;
                           re4    = vc[j+12].real();
                           xre[j+12] = re4;
                           im4    = vc[j+12].imag();
                           xim[j+12] = im4;
                           re5    = vc[j+13].real();
                           xre[j+13] = re5;
                           im5    = vc[j+13].imag();
                           xim[j+13] = im5;
                           re6    = vc[j+14].real();
                           xre[j+14] = re6;
                           im6    = vc[j+14].imag();
                           xim[j+14] = im6;
                           re7    = vc[j+15].real();
                           xre[j+15] = re7;
                           im7    = vc[j+15].imag();
                           xim[j+15] = im7;
                      }
                 }







      } // math


} // gms















#endif /*__GMS_COMPLEX_XMM4R4_HPP__*/
