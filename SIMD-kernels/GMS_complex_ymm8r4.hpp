
#ifndef __GMS_COMPLEX_YMM8C4_HPP__
#define __GMS_COMPLEX_YMM8C4_HPP__ 181020231557


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

    const unsigned int GMS_COMPLEX_YMM8C4_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_YMM8C4_MINOR = 0U;
    const unsigned int GMS_COMPLEX_YMM8C4_MICRO = 0U;
    const unsigned int GMS_COMPLEX_YMM8C4_FULLVER =
      1000U*GMS_COMPLEX_YMM8C4_MAJOR+
      100U*GMS_COMPLEX_YMM8C4_MINOR+
      10U*GMS_COMPLEX_YMM8C4_MICRO;
    const char * const GMS_COMPLEX_YMM8C4_CREATION_DATE = "18-10-2023 15:57 AM +00200 (WED 18 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_YMM8C4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_YMM8C4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_YMM8C4_DESCRIPTION   = "AVX/AVX2 optimized complex number implementation."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace  gms {


       namespace math {
       
       
                   struct __ATTR_ALIGN__(32) ymm8c4_t {
                   
                          __m256 re;
                          __m256 im;
                   };
                   
                   

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_ps(&xre[0]);
                        ymm1  = _mm256_loadu_ps(&yre[0]);
                        _mm256_storeu_ps(&zre[0], _mm256_add_ps(ymm0,ymm1));
                        ymm2  = _mm256_loadu_ps(&xim[0]);
                        ymm3  = _mm256_loadu_ps(&yim[0]);
                        _mm256_storeu_ps(&zim[0], _mm256_add_ps(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       const float * __restrict __ATTR_ALIGN__(32) yre,
                                       const float * __restrict __ATTR_ALIGN__(32) yim,
                                       float *       __restrict __ATTR_ALIGN__(32) zre,
                                       float *       __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_ps(&xre[0]);
                        ymm1  = _mm256_load_ps(&yre[0]);
                        _mm256_store_ps(&zre[0], _mm256_add_ps(ymm0,ymm1));
                        ymm2  = _mm256_load_ps(&xim[0]);
                        ymm3  = _mm256_load_ps(&yim[0]);
                        _mm256_store_ps(&zim[0], _mm256_add_ps(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 yre,
                                     const __m256 yim,
                                     __m256 * __restrict zre,
                                     __m256 * __restrict zim) {
                     
                        register __m256 ymm0,ymm1;
                        ymm0  = _mm256_add_ps(xre,yre);
                        *zre  = ymm0;
                        ymm1  = _mm256_add_ps(xim,yim);
                        *zim  = ymm1;
                }
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cadd_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 yre,
                                          const __m256 yim) {
                                     
                        ymm8c4_t cv;
                        register __m256 ymm0,ymm1;
                        ymm0   = _mm256_add_ps(xre,yre);
                        cv.re  = ymm0;
                        ymm1   = _mm256_add_ps(xim,yim);
                        cv.im  = ymm1;  
                        return (cv);            
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cadd_ymm8c4(const ymm8c4_t x,
                                          const ymm8c4_t y) {
                                     
                        ymm8c4_t cv;
                        register __m256 ymm0,ymm1;
                        ymm0   = _mm256_add_ps(x.re,y.re);
                        cv.re  = ymm0;
                        ymm1   = _mm256_add_ps(x.im,y.im);
                        cv.im  = ymm1;  
                        return (cv);            
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 s,
                                     __m256 * __restrict     zre,
                                     __m256 * __restrict     zim) {

                        *zre = _mm256_add_ps(xre,s);
                        *zim = _mm256_setzero_ps();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ymm8c4_t cadd_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 s) {
                      
                      ymm8c4_t cv;
                      cv.re =  _mm256_add_ps(xre,s);
                      cv.im =  _mm256_setzero_ps();
                      return (cv);                       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ymm8c4_t cadd_ymm8c4(const ymm8c4_t x,
                                          const __m256 s) {
                      
                      ymm8c4_t cv;
                      cv.re =  _mm256_add_ps(x.re,s);
                      cv.im =  _mm256_setzero_ps();
                      return (cv);                       
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm8c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_ps(&xre[0]);
                        ymm1  = _mm256_loadu_ps(&xim[0]);
                        ymm2  = _mm256_loadu_ps(&zre[0]);
                        ymm3  = _mm256_loadu_ps(&zim[0])
                        _mm256_storeu_ps(&zre[0], _mm256_add_ps(ymm2,ymm0));
                        _mm256_storeu_ps(&zim[0], _mm256_add_ps(ymm3,ymm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm8c4_aip(const float * __restrict __ATTR_ALIGN__(32) xre,
                                         const float * __restrict __ATTR_ALIGN__(32) xim,
                                         float *       __restrict __ATTR_ALIGN__(32) zre,
                                         float *       __restrict __ATTR_ALIGN__(32) zim) {
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_ps(&xre[0]);
                        ymm1  = _mm256_load_ps(&xim[0]);
                        ymm2  = _mm256_load_ps(&zre[0]);
                        ymm3  = _mm256_load_ps(&zim[0])
                        _mm256_store_ps(&zre[0], _mm256_add_ps(ymm2,ymm0));
                        _mm256_store_ps(&zim[0], _mm256_add_ps(ymm3,ymm1));
              }


                ////////////////////////////////////////////////////////////////////

                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_ps(&xre[0]);
                        ymm1  = _mm256_loadu_ps(&yre[0]);
                        _mm256_storeu_ps(&zre[0], _mm256_sub_ps(ymm0,ymm1));
                        ymm2  = _mm256_loadu_ps(&xim[0]);
                        ymm3  = _mm256_loadu_ps(&yim[0]);
                        _mm256_storeu_ps(&zim[0], _mm256_sub_ps(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm8c4_a( const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       const float * __restrict __ATTR_ALIGN__(32) yre,
                                       const float * __restrict __ATTR_ALIGN__(32) yim,
                                       float *       __restrict __ATTR_ALIGN__(32) zre,
                                       float *       __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_ps(&xre[0]);
                        ymm1  = _mm256_load_ps(&yre[0]);
                        _mm256_store_ps(&zre[0], _mm256_sub_ps(ymm0,ymm1));
                        ymm2  = _mm256_load_ps(&xim[0]);
                        ymm3  = _mm256_load_ps(&yim[0]);
                        _mm256_store_ps(&zim[0], _mm256_sub_ps(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 yre,
                                     const __m256 yim,
                                     __m256 * __restrict     zre,
                                     __m256 * __restrict     zim) {
                     
                        register __m256 ymm0,ymm1;
                        ymm0  = _mm256_sub_ps(xre,yre);
                        *zre  = ymm0;
                        ymm1  = _mm256_sub_ps(xim,yim);
                        *zim  = ymm1;
                }
                
                
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t csub_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 yre,
                                          const __m256 yim) {
                                    
                        ymm8c4_t cv;
                        register __m256 ymm0,ymm1;
                        ymm0  = _mm256_sub_ps(xre,yre);
                        cv.re  = ymm0;
                        ymm1  = _mm256_sub_ps(xim,yim);
                        cv.im  = ymm1;
                        return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t csub_ymm8c4(const ymm8c4_t x,
                                          const ymm8c4_t y) {
                                    
                        ymm8c4_t cv;
                        register __m256 ymm0,ymm1;
                        ymm0  = _mm256_sub_ps(x.re,y.re);
                        cv.re  = ymm0;
                        ymm1  = _mm256_sub_ps(x.im,y.im);
                        cv.im  = ymm1;
                        return (cv);
                }
                
                


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 s,
                                     __m256 * __restrict     zre,
                                     __m256 * __restrict     zim) {

                        *zre = _mm256_sub_ps(xre,s);
                        *zim = _mm256_setzero_ps();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t csub_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 s) {
                                    
                        ymm8c4_t cv;
                        cv.re = _mm256_sub_ps(xre,s);
                        cv.im = _mm256_setzero_ps();
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t csub_ymm8c4(const ymm8c4_t x,
                                          const __m256 s) {
                                    
                        ymm8c4_t cv;
                        cv.re = _mm256_sub_ps(x.re,s);
                        cv.im = _mm256_setzero_ps();
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm8c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_ps(&xre[0]);
                        ymm1  = _mm256_loadu_ps(&xim[0]);
                        ymm2  = _mm256_loadu_ps(&zre[0]);
                        ymm3  = _mm256_loadu_ps(&zim[0])
                        _mm256_storeu_ps(&zre[0], _mm256_sub_ps(ymm2,ymm0));
                        _mm256_storeu_ps(&zim[0], _mm256_sub_ps(ymm3,ymm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm8c4_aip(const float * __restrict __ATTR_ALIGN__(32) xre,
                                         const float * __restrict __ATTR_ALIGN__(32) xim,
                                         float *       __restrict __ATTR_ALIGN__(32) zre,
                                         float *       __restrict __ATTR_ALIGN__(32) zim) {
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_ps(&xre[0]);
                        ymm1  = _mm256_load_ps(&xim[0]);
                        ymm2  = _mm256_load_ps(&zre[0]);
                        ymm3  = _mm256_load_ps(&zim[0])
                        _mm256_store_ps(&zre[0], _mm256_sub_ps(ymm2,ymm0));
                        _mm256_store_ps(&zim[0], _mm256_sub_ps(ymm3,ymm1));
              }


               ////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                           register __m256 ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_loadu_ps(&xre[0]);
                           ymm1  = _mm256_loadu_ps(&yre[0]);
                           ymm2  = _mm256_loadu_ps(&xim[0]);
                           ymm3  = _mm256_loadu_ps(&yim[0]);
                           ymm4  = _mm256_sub_ps(_mm256_mul_ps(ymm0,ymm1),
                                                                        _mm256_mul_ps(ymm2,ymm3));
                           _mm256_storeu_ps(&zre[0], ymm4);
                           ymm5  = _mm256_mul_ps(_mm256_mul_ps(ymm2,ymm1),
                                                                        _mm256_mul_ps(ymm0,ymm3));
                           _mm256_storeu_ps(&zim[0], ymm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       const float * __restrict __ATTR_ALIGN__(32) yre,
                                       const float * __restrict __ATTR_ALIGN__(32) yim,
                                       float *       __restrict __ATTR_ALIGN__(32) zre,
                                       float *       __restrict __ATTR_ALIGN__(32) zim) {

                           register __m256 ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_load_ps(&xre[0]);
                           ymm1  = _mm256_load_ps(&yre[0]);
                           ymm2  = _mm256_load_ps(&xim[0]);
                           ymm3  = _mm256_load_ps(&yim[0]);
                           ymm4  = _mm256_sub_ps(_mm256_mul_ps(ymm0,ymm1),
                                                                        _mm256_mul_ps(ymm2,ymm3));
                           _mm256_store_ps(&zre[0], ymm4);
                           ymm5  = _mm256_mul_ps(_mm256_mul_ps(ymm2,ymm1),
                                                                        _mm256_mul_ps(ymm0,ymm3));
                           _mm256_store_ps(&zim[0], ymm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 yre,
                                     const __m256 yim,
                                     __m256 * __restrict     zre,
                                     __m256 * __restrict     zim) {

                         register __m256 ymm0,ymm1;
                         ymm0 = _mm256_sub_ps(_mm256_mul_ps(xre,yre),
                                              _mm256_mul_ps(xim,yim));
                         *zre  = ymm0;
                         ymm1 = _mm256_mul_ps(_mm256_mul_ps(xim,yre),
                                              _mm256_mul_ps(xre,yim));
                         *zim  = ymm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cmul_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 yre,
                                          const __m256 yim) {
                                     
                         ymm8c4_t cv
                         register __m256 ymm0,ymm1;
                         ymm0 = _mm256_sub_ps(_mm256_mul_ps(xre,yre),
                                              _mm256_mul_ps(xim,yim));
                         cv.re  = ymm0;
                         ymm1 = _mm256_mul_ps(_mm256_mul_ps(xim,yre),
                                              _mm256_mul_ps(xre,yim));
                         cv.im  = ymm1;
                         return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cmul_ymm8c4(const ymm8c4_t x,
                                          const ymm8c4_t y) {
                                     
                         ymm8c4_t cv
                         register __m256 ymm0,ymm1;
                         ymm0 = _mm256_sub_ps(_mm256_mul_ps(x.re,y.re),
                                              _mm256_mul_ps(x.im,y.im));
                         cv.re  = ymm0;
                         ymm1 = _mm256_mul_ps(_mm256_mul_ps(x.im,y.re),
                                              _mm256_mul_ps(x.re,y.im));
                         cv.im  = ymm1;
                         return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 s,
                                     __m256 * __restrict   zre,
                                     __m256 * __restrict   zim) {

                        *zre = _mm256_mul_ps(xre,s);
                        *zim = _mm256_mul_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cmul_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 s) {
                                     
                        ymm8c4_t cv;
                        cv.re = _mm256_mul_ps(xre,s);
                        cv.im = _mm256_mul_ps(xim,s);
                        return (cv);
               }
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cmul_ymm8c4(const ymm8c4_t x,
                                          const __m256 s) {
                                     
                        ymm8c4_t cv;
                        cv.re = _mm256_mul_ps(x.re,s);
                        cv.im = _mm256_mul_ps(x.im,s);
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm8c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {

                           register __m256 ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_loadu_ps(&xre[0]);
                           ymm1  = _mm256_loadu_ps(&zre[0]);
                           ymm2  = _mm256_loadu_ps(&xim[0]);
                           ymm3  = _mm256_loadu_ps(&zim[0]);
                           ymm4  = _mm256_sub_ps(_mm256_mul_ps(ymm0,ymm1),
                                                 _mm256_mul_ps(ymm2,ymm3));
                           _mm256_storeu_ps(&zre[0], ymm4);
                           ymm5  = _mm256_mul_ps(_mm256_mul_ps(ymm2,ymm1),
                                                 _mm256_mul_ps(ymm0,ymm3));
                           _mm256_storeu_ps(&zim[0], ymm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm8c4_aip(const float * __restrict __ATTR_ALIGN__(32) xre,
                                         const float * __restrict __ATTR_ALIGN__(32) xim,
                                         float *       __restrict __ATTR_ALIGN__(32) zre,
                                         float *       __restrict __ATTR_ALIGN__(32) zim) {

                           register __m256 ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_load_ps(&xre[0]);
                           ymm1  = _mm256_load_ps(&zre[0]);
                           ymm2  = _mm256_load_ps(&xim[0]);
                           ymm3  = _mm256_load_ps(&zim[0]);
                           ymm4  = _mm256_sub_ps(_mm256_mul_ps(ymm0,ymm1),
                                                 _mm256_mul_ps(ymm2,ymm3));
                           _mm256_store_ps(&zre[0], ymm4);
                           ymm5  = _mm256_mul_ps(_mm256_mul_ps(ymm2,ymm1),
                                                 _mm256_mul_ps(ymm0,ymm3));
                           _mm256_store_ps(&zim[0], ymm5);
               }

                 ////////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3; 
                        register __m256 ymm4,ymm5,ymm6;
                        ymm0  = _mm256_loadu_ps(&xre[0]); //a
                        ymm1  = _mm256_loadu_ps(&yim[0]); //d
                        ymm2  = _mm256_loadu_ps(&xim[0]); //b
                        ymm3  = _mm256_loadu_ps(&yre[0]); //c
                        ymm4  = _mm256_fmadd_ps(ymm0,ymm3,
                                                _mm256_mul_ps(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_ps(ymm2,ymm3,
                                                _mm256_mul_ps(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_ps(ymm3,ymm3),
                                                _mm256_mul_ps(ymm1,ymm1));
                        _mm256_storeu_ps(&zre[0], _mm256_div_ps(ymm4,ymm6));
                        _mm256_storeu_ps(&zim[0], _mm256_div_ps(ymm5,ymm6));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4_a(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3; 
                        register __m256 ymm4,ymm5,ymm6;
                        ymm0  = _mm256_load_ps(&xre[0]); //a
                        ymm1  = _mm256_load_ps(&yim[0]); //d
                        ymm2  = _mm256_load_ps(&xim[0]); //b
                        ymm3  = _mm256_load_ps(&yre[0]); //c
                        ymm4  = _mm256_fmadd_ps(ymm0,ymm3,
                                                _mm256_mul_ps(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_ps(ymm2,ymm3,
                                                _mm256_mul_ps(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_ps(ymm3,ymm3,
                                                _mm256_mul_ps(ymm1,ymm1));
                        _mm256_store_ps(&zre[0], _mm256_div_ps(ymm4,ymm6));
                        _mm256_store_ps(&zim[0], _mm256_div_ps(ymm5,ymm6));
                }
              
                                         
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 yre,
                                     const __m256 yim,
                                     __m256 * __restrict zre,
                                     __m256 * __restrict zim) {

                      register __m256 ymm0,ymm1,ymm2;
                      ymm0 = _mm256_fmadd_ps(xre,yre,
                                           _mm256_mul_ps(xim,yim));
                      ymm1 = _mm256_fmsub_ps(xim,yre,
                                           _mm256_mul_ps(xre,yim));
                      ymm2 = _mm256_fmadd_ps(ymm3,ymm3,
                                           _mm256_mul_ps(ymm1,ymm1));
                      *zre  = _mm256_div_ps(ymm0,ymm2);
                      *zim  = _mm256_div_ps(ymm1,ymm2);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 yre,
                                          const __m256 yim) {
                                     
                      ymm8c4_t
                      register __m256 ymm0,ymm1,ymm2;
                      ymm0 = _mm256_fmadd_ps(xre,yre,
                                           _mm256_mul_ps(xim,yim));
                      ymm1 = _mm256_fmsub_ps(x.im,y.re,
                                           _mm256_mul_ps(xre,yim));
                      ymm2 = _mm256_fmadd_ps(ymm3,ymm3,
                                           _mm256_mul_ps(ymm1,ymm1));
                      cv.re  = _mm256_div_ps(ymm0,ymm2);
                      cv.im  = _mm256_div_ps(ymm1,ymm2);
                      return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_ymm8c4(const ymm8c4_t x,
                                          const ymm8c4_t y) {
                                     
                      ymm8c4_t
                      register __m256 ymm0,ymm1,ymm2;
                      ymm0 = _mm256_fmadd_ps(x.re,y.re,
                                           _mm256_mul_ps(x.im,y.im));
                      ymm1 = _mm256_fmsub_ps(x.im,y.re,
                                           _mm256_mul_ps(x.re,y.im));
                      ymm2 = _mm256_fmadd_ps(ymm3,ymm3,
                                           _mm256_mul_ps(ymm1,ymm1));
                      cv.re  = _mm256_div_ps(ymm0,ymm2);
                      cv.im  = _mm256_div_ps(ymm1,ymm2);
                      return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     const __m256 s,
                                     __m256 * __restrict zre,
                                     __m256 * __restrict zim) {

                        *zre = _mm256_div_ps(xre,s);
                        *zim = _mm256_div_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_ymm8c4(const __m256 xre,
                                          const __m256 xim,
                                          const __m256 s) {
                                     
                         ymm8c4_t cv;
                         cv.re = _mm256_div_ps(xre,s);
                         cv.im = _mm256_div_ps(xim,s);
                         return (cv);
               }
               
               
                 
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_ymm8c4(const ymm8c4_t x,
                                          const __m256 s) {
                                     
                         ymm8c4_t cv;
                         cv.re = _mm256_div_ps(x.re,s);
                         cv.im = _mm256_div_ps(x.im,s);
                         return (cv);
               }
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4_s(const __m256 s,
                                       const __m256 xre,
                                       const __m256 xim,
                                       __m256 * __restrict zre,
                                       __m256 * __restrict zim) {
                        
                        register __m256 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm256_setzero_ps();
                        cdiv_ymm8c4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_ymm8c4_s(const __m256 s,
                                            const __m256 xre,
                                            const __m256 xim) {
                                       
                        ymm8c4_t cv;
                        register __m256 t0r,t0i;
                        t0r = s;
                        t0i = _mm256_setzero_ps();
                        cdiv_ymm8c4(t0r,t0i,xre,xim,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_ymm8c4_s(const __m256 s,
                                            const ymm8c4_t x) {
                                       
                        ymm8c4_t cv;
                        register __m256 t0r,t0i;
                        t0r = s;
                        t0i = _mm256_setzero_ps();
                        cdiv_ymm8c4(t0r,t0i,x.re,x.im,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm8c4_s(const __m256 s,
                                             const __m256 xre,
                                             const __m256 xim,
                                             __m256 * __restrict zre,
                                             __m256 * __restrict zim) {
                                             
                        register __m256 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm256_setzero_ps(); 
                        cdiv_smith_ymm8c4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                   
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_smith_ymm8c4_s(const __m256 s,
                                                  const __m256 xre,
                                                  const __m256 xim) {
                                             
                        ymm8c4_t cv;                    
                        register __m256 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm256_setzero_ps(); 
                        cdiv_smith_ymm8c4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_smith_ymm8c4_s(const __m256 s,
                                                  const ymm8c4_t x) {
                                             
                        ymm8c4_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm256_setzero_ps(); 
                        cv = cdiv_smith_ymm8c4(t0,x);
                        return (cv);                 
                 }
                 
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3; 
                        register __m256 ymm4,ymm5,ymm6;
                        ymm0  = _mm256_loadu_ps(&xre[0]); //a
                        ymm1  = _mm256_loadu_ps(&zim[0]); //d
                        ymm2  = _mm256_loadu_ps(&xim[0]); //b
                        ymm3  = _mm256_loadu_ps(&zre[0]); //c
                        ymm4  = _mm256_fmadd_ps(ymm0,ymm3,
                                                _mm256_mul_ps(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_ps(ymm2,ymm3,
                                                _mm256_mul_ps(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_ps(ymm3,ymm3,
                                                _mm256_mul_ps(ymm1,ymm1));
                        _mm256_storeu_ps(&zre[0], _mm256_div_ps(ymm4,ymm6));
                        _mm256_storeu_ps(&zim[0], _mm256_div_ps(ymm5,ymm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm8c4_aip(const float * __restrict __ATTR_ALIGN__(32) xre,
                                         const float * __restrict __ATTR_ALIGN__(32) xim,
                                         float *       __restrict __ATTR_ALIGN__(32) zre,
                                         float *       __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3; 
                        register __m256 ymm4,ymm5,ymm6;
                        ymm0  = _mm256_load_ps(&xre[0]); //a
                        ymm1  = _mm256_load_ps(&zim[0]); //d
                        ymm2  = _mm256_load_ps(&xim[0]); //b
                        ymm3  = _mm256_load_ps(&zre[0]); //c
                        ymm4  = _mm256_fmadd_ps(ymm0,ymm3,
                                                _mm256_mul_ps(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_ps(ymm2,ymm3,
                                                _mm256_mul_ps(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_ps(ymm3,ymm3,
                                                _mm256_mul_ps(ymm1,ymm1));
                        _mm256_store_ps(&zre[0], _mm256_div_ps(ymm4,ymm6));
                        _mm256_store_ps(&zim[0], _mm256_div_ps(ymm5,ymm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm8c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float *       __restrict zre,
                                             float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 r,den;
                        __mmask8 m = 0x0;
                        ymm0 = _mm256_loadu_ps(&yre[0]); // c
                        ymm1 = _mm256_loadu_ps(&yim[0]); // d
                        ymm2 = _mm256_loadu_ps(&xre[0]); // a
                        ymm3 = _mm256_loadu_ps(&xim[0]); // b
                        m    = _mm256_cmp_ps_mask(_mm256_abs_ps(ymm0),
                                                  _mm256_abs_ps(ymm1),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_ps(m,_mm256_div_ps(ymm0,ymm1),
                                                      _mm256_div_ps(ymm1,ymm0)); // r
                        den  = _mm256_mask_blend_ps(m,_mm256_fmadd_ps(r,ymm0,ymm1),
                                                      _mm256_fmadd_ps(r,ymm1,ymm0));
                        _mm256_storeu_ps(&zre[0], _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmadd_ps(ymm2,r,ymm3),den),
                                                _mm256_div_ps(_mm256_fmadd_ps(ymm3,r,ymm2),den)));
                        _mm256_storeu_ps(&zim[0], _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmsub_ps(ymm3,r,ymm2),den),
                                                _mm256_div_ps(_mm256_sub_ps(ymm3,_mm256_mul_ps(ymm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm8c4_a(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float *       __restrict zre,
                                             float *       __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 r,den;
                        __mmask8 m = 0x0;
                        ymm0 = _mm256_load_ps(&yre[0]); // c
                        ymm1 = _mm256_load_ps(&yim[0]); // d
                        ymm2 = _mm256_load_ps(&xre[0]); // a
                        ymm3 = _mm256_load_ps(&xim[0]); // b
                        m    = _mm256_cmp_ps_mask(_mm256_abs_ps(ymm0),
                                                  _mm256_abs_ps(ymm1),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_ps(m,_mm256_div_ps(ymm0,ymm1),
                                                      _mm256_div_ps(ymm1,ymm0)); // r
                        den  = _mm256_mask_blend_ps(m,_mm256_fmadd_ps(r,ymm0,ymm1),
                                                      _mm256_fmadd_ps(r,ymm1,ymm0));
                        _mm256_storeu_ps(&zre[0], _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmadd_ps(ymm2,r,ymm3),den),
                                                _mm256_div_ps(_mm256_fmadd_ps(ymm3,r,ymm2),den)));
                        _mm256_storeu_ps(&zim[0], _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmsub_ps(ymm3,r,ymm2),den),
                                                _mm256_div_ps(_mm256_sub_ps(ymm3,_mm256_mul_ps(ymm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm8c4(const __m256 xre,
                                           const __m256 xim,
                                           const __m256 yre,
                                           const __m256 yim,
                                           __m256 * __restrict zre,
                                           __m256 * __restrict zim) {

                        register __m256 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_ps_mask(_mm256_abs_ps(yre),
                                                  _mm256_abs_ps(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_ps(m,_mm256_div_ps(yre,yim),
                                                      _mm256_div_ps(yim,yre)); // r
                        den  = _mm256_mask_blend_ps(m,_mm256_fmadd_ps(r,yre,yim),
                                                      _mm256_fmadd_ps(r,yim,yre));
                        *zre  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmadd_ps(xre,r,xim),den),
                                                _mm256_div_ps(_mm256_fmadd_ps(xim,r,xre),den));
                        *zim  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmsub_ps(xim,r,xre),den),
                                                _mm256_div_ps(_mm256_sub_ps(xim,_mm256_mul_ps(xre,r)),den)));
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_smith_ymm8c4(const __m256 xre,
                                                const __m256 xim,
                                                const __m256 yre,
                                                const __m256 yim) {
                                           
                        ymm8c4_t cv
                        register __m256 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_ps_mask(_mm256_abs_ps(yre),
                                                  _mm256_abs_ps(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_ps(m,_mm256_div_ps(yre,yim),
                                                      _mm256_div_ps(yim,yre)); // r
                        den  = _mm256_mask_blend_ps(m,_mm256_fmadd_ps(r,yre,yim),
                                                      _mm256_fmadd_ps(r,yim,yre));
                        cv.re  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmadd_ps(xre,r,xim),den),
                                                _mm256_div_ps(_mm256_fmadd_ps(xim,r,xre),den));
                        cv.im  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmsub_ps(xim,r,xre),den),
                                                _mm256_div_ps(_mm256_sub_ps(xim,_mm256_mul_ps(xre,r)),den)));
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cdiv_smith_ymm8c4(const ymm8c4_t x,
                                                const ymm8c4_t y) {
                                           
                        ymm8c4_t cv
                        register __m256 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_ps_mask(_mm256_abs_ps(y.re),
                                                  _mm256_abs_ps(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_ps(m,_mm256_div_ps(y.re,y.im),
                                                      _mm256_div_ps(y.im,y.re)); // r
                        den  = _mm256_mask_blend_ps(m,_mm256_fmadd_ps(r,y.re,y.im),
                                                      _mm256_fmadd_ps(r,y.im,y.re));
                        cv.re  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmadd_ps(x.re,r,x.im),den),
                                                _mm256_div_ps(_mm256_fmadd_ps(x.im,r,x.re),den));
                        cv.im  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmsub_ps(x.im,r,x.re),den),
                                                _mm256_div_ps(_mm256_sub_ps(x.im,_mm256_mul_ps(x.re,r)),den)));
                        return (cv);
               }





                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_ymm8c4_u(const float * __restrict re,
                                       const float * __restrict im,
                                       float * __restrict  cabs) {

                        register __m256 ymm0,ymm1,ymm2,ymm3,ymm4;
                        ymm0  = _mm256_loadu_ps(&re[0]);
                        ymm1  = _mm256_mul_ps(ymm0,ymm0);
                        ymm2  = _mm256_loadu_ps(&im[0]);
                        ymm3  = _mm256_mul_ps(ymm2,ymm2);
                        ymm4  = _mm256_sqrt_ps(_mm256_add_ps(ymm1,ymm3));
                        _mm256_storeu_ps(&cabs[0],ymm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) re,
                                       const float * __restrict __ATTR_ALIGN__(32) im,
                                       float * __restrict  __ATTR_ALIGN__(32) cabs) {

                        register __m256 ymm0,ymm1,ymm2,ymm3,ymm4;
                        ymm0  = _mm256_load_ps(&re[0]);
                        ymm1  = _mm256_mul_ps(ymm0,ymm0);
                        ymm2  = _mm256_load_ps(&im[0]);
                        ymm3  = _mm256_mul_ps(ymm2,ymm2);
                        ymm4  = _mm256_sqrt_ps(_mm256_add_ps(ymm1,ymm3));
                        _mm256_store_ps(&cabs[0],ymm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256 cabs_ymm8c4(const __m256 re,
                                       const __m256 im) {

                        register __m256 ymm0,ymm1,cabs;
                        ymm0 = _mm256_mul_ps(re,re);
                        ymm1 = _mm256_mul_ps(im,im);
                        cabs = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        return (cabs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256 cabs_ymm8c4(const ymm8c4_t x) {

                        register __m256 ymm0,ymm1,cabs;
                        ymm0 = _mm256_mul_ps(x.re,x.re);
                        ymm1 = _mm256_mul_ps(x.im,x.im);
                        cabs = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        return (cabs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_ymm8c4_u(const float * __restrict re,
                                       const float * __restrict im,
                                       float * __restrict  carg) {

                        register __m256 ymm0,ymm1;
                        ymm0 = _mm256_loadu_ps(&re[0]);
                        ymm1 = _mm256_loadu_ps(&im[0]);
                        _mm256_storeu_ps(&carg[0], xatan2f(ymm0,ymm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) re,
                                       const float * __restrict __ATTR_ALIGN__(32) im,
                                       float * __restrict  __ATTR_ALIGN__(32) carg) {

                        register __m256 ymm0,ymm1;
                        ymm0 = _mm256_load_ps(&re[0]);
                        ymm1 = _mm256_load_ps(&im[0]);
                        _mm256_store_ps(&carg[0], xatan2f(ymm0,ymm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256 carg_ymm8c4(const __m256 re,
                                       const __m256 im) {

                       register __m256 carg;
                       carg = _mm256_atan2_ps(re,im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256 carg_ymm8c4(ymm8c4_t x) {

                       register __m256 carg;
                       carg = _mm256_atan2_ps(x.re,x.im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_ymm8c4(const __m256 re,
	                             const __m256 im,
	                             __m256 * __restrict clogr,
	                             __m256 * __restrict clogi) {
	                
	                register __m256 t1,t2,ln;
	                t1  = cabs_ymm8c4(re,im);
	                t2  = carg_ymm8c4(re,im);
	                ln  = _mm256_log_ps(t1);
	                *clogr = ln;
	                *clogi = t2;                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) pre,
	                               const float * __restrict __ATTR_ALIGN__(32) pim,
	                               float * __restrict clogr,
	                               float * __restrict clogi) {
	                
	                register __m256 re = _mm256_load_ps(&pre[0]);
	                register __m256 im = _mm256_load_ps(&pim[0]);
	                register __m256 t1,t2,ln;
	                t1  = cabs_ymm8c4(re,im);
	                t2  = carg_ymm8c4(re,im);
	                ln  = _mm256_log_ps(t1);
	                _mm256_store_ps(&clogr[0] ,ln);
	                _mm256_store_ps(&clogi[0] ,t2);                    
	        }
	        
	        
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_ymm8c4_u(const float * __restrict  pre,
	                               const float * __restrict  pim,
	                               float * __restrict clogr,
	                               float * __restrict clogi) {
	                
	                register __m256 re = _mm256_loadu_ps(&pre[0]);
	                register __m256 im = _mm256_loadu_ps(&pim[0]);
	                register __m256 t1,t2,ln;
	                t1  = cabs_ymm8c4(re,im);
	                t2  = carg_ymm8c4(re,im);
	                ln  = _mm256_log_ps(t1);
	                _mm256_storeu_ps(&clogr[0] ,ln);
	                _mm256_storeu_ps(&clogi[0] ,t2);                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ymm8c4_t clog_ymm8c4(const ymm8c4_t x){
	                                  
	                ymm8c4_t clog;                           
	                register __m256 t1,t2,ln;
	                t1  = cabs_ymm8c4(x.re,x.im);
	                t2  = carg_ymm8c4(x.re,x.im);
	                ln  = _mm256_log_ps(t1);
	                clog.re = ln;
	                clog.im = t2;
	                return (clog);                    
	        }

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm8c4_u(float * __restrict re,
                                        float * __restrict im) {

                        register __m256 c;
                        c = negate_ymm8r4(_mm256_loadu_ps(&im[0]));
                        _mm256_storeu_ps(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm8c4_a(float * __restrict __ATTR_ALIGN__(32) re,
                                        float * __restrict __ATTR_ALIGN__(32) im) {
                                        
                        register __m256 c;
                        c = negate_ymm8r4(_mm256_load_ps(&im[0]));
                        _mm256_store_ps(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm8c4(__m256 * __restrict re,
                                      __m256 * __restrict im) {
                         
                        register __m256 c;              
                        c = negate_ymm8r4(*im);
                        *im = c;
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm8c4_v2(const __m256 xre,
                                         const __m256 xim,
                                         __m256 * __restrict yre,
                                         __m256 * __restrict yim) {
                         
                        //register __m256 c;              
                        //c = negate_ymm8c4(*im);
                        //*im = c;
                        *yre = xre; 
                        *yim = negate_ymm8r4(xim);
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cconj_ymm8c4_v2(const __m256 xre,
                                              const __m256 xim) {                                              
                         
                        //register __m256 c;              
                        //c = negate_ymm8c4(*im);
                        //*im = c;
                        ymm8c4_t cv;
                        cv.re = xre; 
                        cv.im = negate_ymm8r4(xim);
                        return (cv);
                   } 
                   
                   
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cconj_ymm8c4_v2(const ymm8c4_t x) {                                              
                         
                        //register __m256 c;              
                        //c = negate_ymm8c4(*im);
                        //*im = c;
                        ymm8c4_t cv;
                        cv.re = x.re; 
                        cv.im = negate_ymm8r4(x.im);
                        return (cv);
                   } 
                   
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict  csre,
                                       float * __restrict  csim) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_loadu_ps(&xre[0]);
                      ymm1  = _mm256_loadu_ps(&xim[0]);
                      ymm2  = _mm256_mul_ps(_mm256_cos_ps(ymm0),_mm256_cosh_ps(ymm1));
                      _mm256_storeu_ps(&csre[0],ymm2);
                      ymm3  = _mm256_mul_ps(_mm256_sin_ps(ymm0),_mm256_sinh_ps(ymm1));
                      _mm256_storeu_ps(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       float * __restrict  __ATTR_ALIGN__(32) csre,
                                       float * __restrict  __ATTR_ALIGN__(32) csim) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_load_ps(&xre[0]);
                      ymm1  = _mm256_load_ps(&xim[0]);
                      ymm2  = _mm256_mul_ps(_mm256_cos_ps(ymm0),_mm256_cosh_ps(ymm1));
                      _mm256_store_ps(&csre[0],ymm2);
                      ymm3  = _mm256_mul_ps(_mm256_sin_ps(ymm0),_mm256_sinh_ps(ymm1));
                      _mm256_store_ps(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     __m256 * __restrict csre,
                                     __m256 * __restrict csim) {

                      register __m256 ymm0,ymm1;
                      ymm0  = _mm256_mul_ps(_mm256_cos_ps(xre),_mm256_cosh_ps(xim));
                      *csre = ymm0;
                      ymm1  = _mm256_mul_ps(_mm256_sin_ps(xre),_mm256_sinh_ps(xim));
                      *csim = ymm1; 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t ccos_ymm8c4(const __m256 xre,
                                          const __m256 xim) {
                                    
                      ymm8c4_t cv;
                      register __m256 ymm0,ymm1;
                      ymm0  = _mm256_mul_ps(_mm256_cos_ps(xre),_mm256_cosh_ps(xim));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_ps(_mm256_sin_ps(xre),_mm256_sinh_ps(xim));
                      cv.im = ymm1;
                      return (cv); 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t ccos_ymm8c4(const ymm8c4_t x) {
                                    
                      ymm8c4_t cv;
                      register __m256 ymm0,ymm1;
                      ymm0  = _mm256_mul_ps(_mm256_cos_ps(x.re),_mm256_cosh_ps(x.im));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_ps(_mm256_sin_ps(x.re),_mm256_sinh_ps(x.im));
                      cv.im = ymm1;
                      return (cv); 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict  csre,
                                       float * __restrict  csim) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_loadu_ps(&xre[0]);
                      ymm1  = _mm256_loadu_ps(&xim[0]);
                      ymm2  = _mm256_mul_ps(_mm256_cosh_ps(ymm0),_mm256_cos_ps(ymm1));
                      _mm256_storeu_ps(&csre[0],ymm2);
                      ymm3  = _mm256_mul_ps(_mm256_sinh_ps(ymm0),_mm256_sin_ps(ymm1));
                      _mm256_storeu_ps(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       float * __restrict  __ATTR_ALIGN__(32) csre,
                                       float * __restrict  __ATTR_ALIGN__(32) csim) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_load_ps(&xre[0]);
                      ymm1  = _mm256_load_ps(&xim[0]);
                      ymm2  = _mm256_mul_ps(_mm256_cosh_ps(ymm0),_mm256_cos_ps(ymm1));
                      _mm256_store_ps(&csre[0],ymm2);
                      ymm3  = _mm256_mul_ps(_mm256_sinh_ps(ymm0),_mm256_sin_ps(ymm1));
                      _mm256_store_ps(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     __m256 * __restrict csre,
                                     __m256 * __restrict csim) {

                      register __m256 ymm0,ymm1;
                      ymm0  = _mm256_mul_ps(_mm256_cosh_ps(xre),_mm256_cos_ps(xim));
                      *csre = ymm0;
                      ymm1  = _mm256_mul_ps(_mm256_sinh_ps(xre),_mm256_sin_ps(xim));
                      *csim = ymm1; 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t ccosh_ymm8c4(const __m256 xre,
                                           const __m256 xim) {
                                          
                      ymm8c4_t cv;
                      register __m256 ymm0,ymm1;
                      ymm0  = _mm256_mul_ps(_mm256_cosh_ps(xre),_mm256_cos_ps(xim));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_ps(_mm256_sinh_ps(xre),_mm256_sin_ps(xim));
                      cv.im = ymm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t ccosh_ymm8c4(const ymm8c4_t x) {
                                          
                      ymm8c4_t cv;
                      register __m256 ymm0,ymm1;
                      ymm0  = _mm256_mul_ps(_mm256_cosh_ps(x.re),_mm256_cos_ps(x.im));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_ps(_mm256_sinh_ps(x.re),_mm256_sin_ps(x.im));
                      cv.im = ymm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_ymm8c4(const __m256 xre,
	                             const __m256 xim,
	                             const float n,
	                             __m256 * __restrict powr,
	                             __m256 * __restrict powi) {
	                             
	                register __m256 ymm0,ymm1;
	                register __m256 r,tht;
	                register __m256 vn,pt;
	                register __m256 ta;
	                ymm0  = _mm256_mul_ps(xre,xre);
	                vn    = _mm256_set1_ps(n);
	                ymm1  = _mm256_mul_ps(xim,xim);
	                r     = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
	                tht   = _mm256_atan_ps(_mm256_div_ps(xim,xre));
	                pt    = _mm256_pow_ps(r,vn);
	                ta    = _mm256_mul_ps(vn,tht);
	                *powr = _mm256_mul_ps(pt,_mm256_cos_ps(ta));
	                *powi = _mm256_mul_ps(pt,_mm256_sin_ps(ta));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) pxre,
	                               const float * __restrict __ATTR_ALIGN__(32) pxim,
	                               const float n,
	                               float * __restrict __ATTR_ALIGN__(32) ppowr
	                               float * __restrict __ATTR_ALIGN__(32) ppowi) {
	                  
	                register __m256 xre = _mm256_load_ps(&pxre[0]);
	                register __m256 xim = _mm256_load_ps(&pxim[0]);           
	                register __m256 ymm0,ymm1;
	                register __m256 r,tht;
	                register __m256 vn,pt;
	                register __m256 ta;
	                ymm0  = _mm256_mul_ps(xre,xre);
	                vn    = _mm256_set1_ps(n);
	                ymm1  = _mm256_mul_ps(xim,xim);
	                r     = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
	                tht   = _mm256_atan_ps(_mm256_div_ps(xim,xre));
	                pt    = _mm256_pow_ps(r,vn);
	                ta    = _mm256_mul_ps(vn,tht);
	                _mm256_store_ps(&ppowr[0] ,_mm256_mul_ps(pt,_mm256_cos_ps(ta)));
	                _mm256_store_ps(&ppowi[0] ,_mm256_mul_ps(pt,_mm256_sin_ps(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_ymm8c4_u(const float * __restrict  pxre,
	                               const float * __restrict  pxim,
	                               const float n,
	                               float * __restrict  ppowr
	                               float * __restrict  ppowi) {
	                  
	                register __m256 xre = _mm256_loadu_ps(&pxre[0]);
	                register __m256 xim = _mm256_loadu_ps(&pxim[0]);           
	                register __m256 ymm0,ymm1;
	                register __m256 r,tht;
	                register __m256 vn,pt;
	                register __m256 ta;
	                ymm0  = _mm256_mul_ps(xre,xre);
	                vn    = _mm256_set1_ps(n);
	                ymm1  = _mm256_mul_ps(xim,xim);
	                r     = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
	                tht   = _mm256_atan_ps(_mm256_div_ps(xim,xre));
	                pt    = _mm256_pow_ps(r,vn);
	                ta    = _mm256_mul_ps(vn,tht);
	                _mm256_storeu_ps(&ppowr[0] ,_mm256_mul_ps(pt,_mm256_cos_ps(ta)));
	                _mm256_storeu_ps(&ppowi[0] ,_mm256_mul_ps(pt,_mm256_sin_ps(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ZMM16c4_t cpow_ymm8c4(const ZMM16c4_t x,
	                                  const float n) {
	                   
	                ZMM16c4_t cp;        
	                register __m256 ymm0,ymm1;
	                register __m256 r,tht;
	                register __m256 vn,pt;
	                register __m256 ta;
	                ymm0  = _mm256_mul_ps(x.re,x.re);
	                vn    = _mm256_set1_ps(n);
	                ymm1  = _mm256_mul_ps(x.im,x.im);
	                r     = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
	                tht   = _mm256_atan_ps(_mm256_div_ps(x.im,x.re));
	                pt    = _mm256_pow_ps(r,vn);
	                ta    = _mm256_mul_ps(vn,tht);
	                cp.re = _mm256_mul_ps(pt,_mm256_cos_ps(ta));
	                cp.im = _mm256_mul_ps(pt,_mm256_sin_ps(ta));      
	                return (cp);              
	       }
	       


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_ymm8c4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask8 * __restrict eqr,
                                      __mmask8 * __restrict eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_loadu_ps(&xre[0]);
                      ymm1 = _mm256_loadu_ps(&yre[0]);
                      _mm256_storeu_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_EQ_OQ));
                      ymm2 = _mm256_loadu_ps(&xim[0]);
                      ymm3 = _mm256_loadu_ps(&yim[0]);
                      _mm256_storeu_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_EQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                      const float * __restrict __ATTR_ALIGN__(32) xim,
                                      const float * __restrict __ATTR_ALIGN__(32) yre,
                                      const float * __restrict __ATTR_ALIGN__(32) yim,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqr,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_load_ps(&xre[0]);
                      ymm1 = _mm256_load_ps(&yre[0]);
                      _mm256_store_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_EQ_OQ));
                      ymm2 = _mm256_load_ps(&xim[0]);
                      ymm3 = _mm256_load_ps(&yim[0]);
                      _mm256_store_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_EQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_ymm8c4(const __m256 xre,
                                    const __m256 xim,
                                    const __m256 yre,
                                    const __m256 yim,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(xre,yre,_CMP_EQ_OQ);
                         *eqi = _mm256_cmp_ps_mask(xim,yim,_CMP_EQ_OQ);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_ymm8c4(const ymm8c4_t x,
                                    const ymm8c4_t y,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_EQ_OQ);
                         *eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_EQ_OQ);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_ymm8c4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask8 * __restrict eqr,
                                      __mmask8 * __restrict eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_loadu_ps(&xre[0]);
                      ymm1 = _mm256_loadu_ps(&yre[0]);
                      _mm256_storeu_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_GT_OQ));
                      ymm2 = _mm256_loadu_ps(&xim[0]);
                      ymm3 = _mm256_loadu_ps(&yim[0]);
                      _mm256_storeu_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_GT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                      const float * __restrict __ATTR_ALIGN__(32) xim,
                                      const float * __restrict __ATTR_ALIGN__(32) yre,
                                      const float * __restrict __ATTR_ALIGN__(32) yim,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqr,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_load_ps(&xre[0]);
                      ymm1 = _mm256_load_ps(&yre[0]);
                      _mm256_store_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_GT_OQ));
                      ymm2 = _mm256_load_ps(&xim[0]);
                      ymm3 = _mm256_load_ps(&yim[0]);
                      _mm256_store_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_GT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_ymm8c4(const __m256 xre,
                                    const __m256 xim,
                                    const __m256 yre,
                                    const __m256 yim,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(xre,yre,_CMP_GT_OQ);
                         *eqi = _mm256_cmp_ps_mask(xim,yim,_CMP_GT_OQ);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_ymm8c4(const ymm8c4_t x,
                                    const ymm8c4_t y,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_GT_OQ);
                         *eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_GT_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_ymm8c4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask8 * __restrict eqr,
                                      __mmask8 * __restrict eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_loadu_ps(&xre[0]);
                      ymm1 = _mm256_loadu_ps(&yre[0]);
                      _mm256_storeu_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_LT_OQ));
                      ymm2 = _mm256_loadu_ps(&xim[0]);
                      ymm3 = _mm256_loadu_ps(&yim[0]);
                      _mm256_storeu_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_LT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                      const float * __restrict __ATTR_ALIGN__(32) xim,
                                      const float * __restrict __ATTR_ALIGN__(32) yre,
                                      const float * __restrict __ATTR_ALIGN__(32) yim,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqr,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_load_ps(&xre[0]);
                      ymm1 = _mm256_load_ps(&yre[0]);
                      _mm256_store_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_LT_OQ));
                      ymm2 = _mm256_load_ps(&xim[0]);
                      ymm3 = _mm256_load_ps(&yim[0]);
                      _mm256_store_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_LT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_ymm8c4(const __m256 xre,
                                    const __m256 xim,
                                    const __m256 yre,
                                    const __m256 yim,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(xre,yre,_CMP_LT_OQ);
                         *eqi = _mm256_cmp_ps_mask(xim,yim,_CMP_LT_OQ);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_ymm8c4(const ymm8c4_t x,
                                    const ymm8c4_t y,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_LT_OQ);
                         *eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_LT_OQ);
              }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_ymm8c4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask8 * __restrict eqr,
                                      __mmask8 * __restrict eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_loadu_ps(&xre[0]);
                      ymm1 = _mm256_loadu_ps(&yre[0]);
                      _mm256_storeu_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_NEQ_OQ));
                      ymm2 = _mm256_loadu_ps(&xim[0]);
                      ymm3 = _mm256_loadu_ps(&yim[0]);
                      _mm256_storeu_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_NEQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                      const float * __restrict __ATTR_ALIGN__(32) xim,
                                      const float * __restrict __ATTR_ALIGN__(32) yre,
                                      const float * __restrict __ATTR_ALIGN__(32) yim,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqr,
                                      __mmask8 * __restrict __ATTR_ALIGN__(32) eqi ) {

                      register __m256 ymm0,ymm1,ymm2,ymm3;
                      ymm0 = _mm256_load_ps(&xre[0]);
                      ymm1 = _mm256_load_ps(&yre[0]);
                      _mm256_store_ps(&eqr[0],
                                       _mm256_cmp_ps_mask(ymm0,ymm1,_CMP_NEQ_OQ));
                      ymm2 = _mm256_load_ps(&xim[0]);
                      ymm3 = _mm256_load_ps(&yim[0]);
                      _mm256_store_ps(&eqi[0],
                                       _mm256_cmp_ps_mask(ymm2,ymm3,_CMP_NEQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_ymm8c4(const __m256 xre,
                                    const __m256 xim,
                                    const __m256 yre,
                                    const __m256 yim,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(xre,yre,_CMP_NEQ_OQ);
                         *eqi = _mm256_cmp_ps_mask(xim,yim,_CMP_NEQ_OQ);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_ymm8c4(const ymm8c4_t x,
                                     const ymm8c4_t y,
                                    __mmask8 * __restrict eqr,
                                    __mmask8 * __restrict eqi) {

                         *eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_NEQ_OQ);
                         *eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_NEQ_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict cexpr,
                                       float * __restrict cexpi ) {

                        register const __m256 I = _mm256_set1_ps(1.0f);
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_ps(&xre[0]);
                        ymm1  = _mm256_loadu_ps(&xim[0]);
                        ymm2  = _mm256_exp_ps(ymm0);
                        ymm3  = _mm256_mul_ps(ymm2,_mm256_cos_ps(ymm1));
                        _mm256_storeu_ps(&cexpr[0],ymm3);
                        ymm4  = _mm256_mul_ps(ymm2,_mm256_mul_ps(_mm256_sin_ps(ymm1),I));
                        _mm256_storeu_ps(&cexpi[0],ymm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       float * __restrict __ATTR_ALIGN__(32) cexpr,
                                       float * __restrict __ATTR_ALIGN__(32) cexpi ) {

                        register const __m256 I = _mm256_set1_ps(1.0f);
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_ps(&xre[0]);
                        ymm1  = _mm256_load_ps(&xim[0]);
                        ymm2  = _mm256_exp_ps(ymm0);
                        ymm3  = _mm256_mul_ps(ymm2,_mm256_cos_ps(ymm1));
                        _mm256_store_ps(&cexpr[0],ymm3);
                        ymm4  = _mm256_mul_ps(ymm2,_mm256_mul_ps(_mm256_sin_ps(ymm1),I));
                        _mm256_store_ps(&cexpi[0],ymm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_ymm8c4(const __m256 xre,
                                     const __m256 xim,
                                     __m256 * __restrict cexpr,
                                     __m256 * __restrict cexpi) {

                        register const __m256 I = _mm256_set1_ps(1.0f);
                        register __m256 ymm0;
                        ymm0   = _mm256_exp_ps(xre);
                        *cexpr = _mm256_mul_ps(ymm0,_mm256_cos_ps(xim));
                        *cexpi = _mm256_mul_ps(ymm0,_mm256_mul_ps(_mm256_sin_ps(xim),I));
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cexp_ymm8c4(const __m256 xre,
                                          const __m256 xim) {
                                     
                        ymm8c4_t cv;
                        register const __m256 I = _mm256_set1_ps(1.0f);
                        register __m256 ymm0;
                        ymm0   = _mm256_exp_ps(xre);
                        cv.re = _mm256_mul_ps(ymm0,_mm256_cos_ps(xim));
                        cv.im = _mm256_mul_ps(ymm0,_mm256_mul_ps(_mm256_sin_ps(xim),I));
                        return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cexp_ymm8c4(const ymm8c4_t x) {
                                     
                        ymm8c4_t cv;
                        register const __m256 I = _mm256_set1_ps(1.0f);
                        register __m256 ymm0;
                        ymm0   = _mm256_exp_ps(x.re);
                        cv.re = _mm256_mul_ps(ymm0,_mm256_cos_ps(x.im));
                        cv.im = _mm256_mul_ps(ymm0,_mm256_mul_ps(_mm256_sin_ps(x.im),I));
                        return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_ymm8c4_u(const float * __restrict rho,
                                         const float * __restrict tht,
                                         float * __restrict  re,
                                         float * __restrict  im) {

                         register __m256 ymm0,ymm1,ymm2,ymm3;
                         ymm0 = _mm256_loadu_ps(&rho[0]);
                         ymm1 = _mm256_loadu_ps(&tht[0]);
                         ymm2 = _mm256_mul_ps(ymm0,_mm256_cos_ps(ymm1)); //tht
                         _mm256_storeu_ps(&re[0],ymm2);
                         ymm3 = _mm256_mul_ps(ymm0,_mm256_sin_ps(ymm1)); //tht
                         _mm256_storeu_ps(&im[0],ymm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) rho,
                                         const float * __restrict __ATTR_ALIGN__(32) tht,
                                         float * __restrict  __ATTR_ALIGN__(32) re,
                                         float * __restrict  __ATTR_ALIGN__(32) im) {

                         register __m256 ymm0,ymm1,ymm2,ymm3;
                         ymm0 = _mm256_load_ps(&rho[0]);
                         ymm1 = _mm256_load_ps(&tht[0]);
                         ymm2 = _mm256_mul_ps(ymm0,_mm256_cos_ps(ymm1)); //tht
                         _mm256_store_ps(&re[0],ymm2);
                         ymm3 = _mm256_mul_ps(ymm0,_mm256_sin_ps(ymm1)); //tht
                         _mm256_store_ps(&im[0],ymm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_ymm8c4(const __m256 rho,
                                       const __m256 tht,
                                       __m256 * __restrict re,
                                       __m256 * __restrict im) {

                        register __m256 ymm0,ymm1;
                        ymm0 = _mm256_mul_ps(rho,_mm256_cos_ps(tht));
                        *re  = ymm0;
                        ymm1 = _mm256_mul_ps(rho,_mm256_sin_ps(tht));
                        *im  = ymm1;
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cpolar_ymm8c4(const __m256 rho,
                                            const __m256 tht) {
                                      
                        ymm8c4_t cv
                        register __m256 ymm0,ymm1;
                        ymm0 = _mm256_mul_ps(rho,_mm256_cos_ps(tht));
                        cv.re  = ymm0;
                        ymm1 = _mm256_mul_ps(rho,_mm256_sin_ps(tht));
                        cv.im  = ymm1;
                        return (cv);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_ymm8c4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict wrkc,
                                       float * __restrict csqr,
                                       float * __restrict csqi) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 half = _mm256_set1_ps(0.5f);
                        cabs_ymm8c4_u(xre,xim,wrkc);
                        ymm0  = _mm256_loadu_ps(&xre[0]);
                        ymm1  = _mm256_loadu_ps(&wrkc[0]);
                        ymm2  = _mm256_mul_ps(half,_mm256_add_ps(ymm1,ymm0));
                        _mm256_storeu_ps(&csqr[0],_mm256_sqrt_ps(ymm2));
                        ymm3  = _mm256_mul_ps(half,_mm256_sub_ps(ymm1,ymm0));
                        _mm256_storeu_ps(&csqi[0],_mm256_sqrt_ps(ymm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                       const float * __restrict __ATTR_ALIGN__(32) xim,
                                       float * __restrict __ATTR_ALIGN__(32) wrkc,
                                       float * __restrict __ATTR_ALIGN__(32) csqr,
                                       float * __restrict __ATTR_ALIGN__(32) csqi) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 half = _mm256_set1_ps(0.5f);
                        cabs_ymm8c4_a(xre,xim,wrkc);
                        ymm0  = _mm256_load_ps(&xre[0]);
                        ymm1  = _mm256_load_ps(&wrkc[0]);
                        ymm2  = _mm256_mul_ps(half,_mm256_add_ps(ymm1,ymm0));
                        _mm256_store_ps(&csqr[0],_mm256_sqrt_ps(ymm2));
                        ymm3  = _mm256_mul_ps(half,_mm256_sub_ps(ymm1,ymm0));
                        _mm256_store_ps(&csqi[0],_mm256_sqrt_ps(ymm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_ymm8c4(const __m256 xre,
                                      const __m256 xim,
                                      __m256 * __restrict wrkc,
                                      __m256 * __restrict csqr,
                                      __m256 * __restrict csqi) {

                       register __m256 ymm0,ymm1;
                       register __m256 half = _mm256_set1_ps(0.5f); 
                       cabs_ymm8c4(xre,xim,wrkc);
                       ymm0  = _mm256_mul_ps(half,_mm256_add_ps(*wrkc,xre));
                       *csqr = ymm0;
                       ymm1  = _mm256_mul_ps(half,_mm256_sub_ps(*wrkc,xre));
                       *csqi = ymm1; 
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t csqrt_ymm8c4(const __m256 xre,
                                           const __m256 xim,
                                          __m256 * __restrict wrkc) {
                                          
                       ymm8c4_t cv;
                       register __m256 ymm0,ymm1;
                       register __m256 half = _mm256_set1_ps(0.5f); 
                       cabs_ymm8c4(xre,xim,wrkc);
                       ymm0  = _mm256_mul_ps(half,_mm256_add_ps(*wrkc,xre));
                       cv.re = ymm0;
                       ymm1  = _mm256_mul_ps(half,_mm256_sub_ps(*wrkc,xre));
                       cv.im = ymm1; 
                       return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t csqrt_ymm8c4(const ymm8c4_t x,
                                          __m256 * __restrict wrkc) {
                                          
                       ymm8c4_t cv;
                       register __m256 ymm0,ymm1;
                       register __m256 half = _mm256_set1_ps(0.5f); 
                       cabs_ymm8c4(x.re,x.im,wrkc);
                       ymm0  = _mm256_mul_ps(half,_mm256_add_ps(*wrkc,x.re));
                       cv.re = ymm0;
                       ymm1  = _mm256_mul_ps(half,_mm256_sub_ps(*wrkc,x.re));
                       cv.im = ymm1; 
                       return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_ymm8c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict zre,
                                             float * __restrict zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        rep  = _mm256_fmsub_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        imp  = _mm256_fmadd_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        ymm0 = _mm256_mul_ps(rep,rep);
                        ymm1 = _mm256_mul_ps(imp,imp);
                        ymm2 = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        _mm256_storeu_ps(&zre[0], _mm256_div_ps(rep,ymm2));
                        _mm256_storeu_ps(&zim[0], _mm256_div_ps(imp,ymm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                             const float * __restrict __ATTR_ALIGN__(32) xim,
                                             const float * __restrict __ATTR_ALIGN__(32) yre,
                                             const float * __restrict __ATTR_ALIGN__(32) yim,
                                             float * __restrict __ATTR_ALIGN__(32) zre,
                                             float * __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        rep  = _mm256_fmsub_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        imp  = _mm256_fmadd_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        ymm0 = _mm256_mul_ps(rep,rep);
                        ymm1 = _mm256_mul_ps(imp,imp);
                        ymm2 = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        _mm256_store_ps(&zre[0], _mm256_div_ps(rep,ymm2));
                        _mm256_store_ps(&zim[0], _mm256_div_ps(imp,ymm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_ymm8c4(  const __m256  xre,
                                             const __m256  xim,
                                             const __m256  yre,
                                             const __m256  yim,
                                             __m256 * __restrict zre,
                                             __m256 * __restrict zim) {

                        register __m256 rep,imp,ymm0,ymm1,ymm2;
                        rep  = _mm256_fmsub_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        imp  = _mm256_fmadd_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        ymm0 = _mm256_mul_ps(rep,rep);
                        ymm1 = _mm256_mul_ps(imp,imp);
                        ymm2 = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        *zre = _mm256_div_ps(rep,ymm2);
                        *zim = _mm256_div_ps(imp,ymm2);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cnorm_prod_ymm8c4(  const __m256  xre,
                                                  const __m256  xim,
                                                  const __m256  yre,
                                                  const __m256  yim) {
                                             
                        ymm8c4_t cv;
                        register __m256 rep,imp,ymm0,ymm1,ymm2;
                        rep  = _mm256_fmsub_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        imp  = _mm256_fmadd_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        ymm0 = _mm256_mul_ps(rep,rep);
                        ymm1 = _mm256_mul_ps(imp,imp);
                        ymm2 = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        cv.re = _mm256_div_ps(rep,ymm2);
                        cv.im = _mm256_div_ps(imp,ymm2);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cnorm_prod_ymm8c4(  const ymm8c4_t x,
                                                  const ymm8c4_t y) {
                                             
                        ymm8c4_t cv;
                        register __m256 rep,imp,ymm0,ymm1,ymm2;
                        rep  = _mm256_fmsub_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        imp  = _mm256_fmadd_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        ymm0 = _mm256_mul_ps(rep,rep);
                        ymm1 = _mm256_mul_ps(imp,imp);
                        ymm2 = _mm256_sqrt_ps(_mm256_add_ps(ymm0,ymm1));
                        cv.re = _mm256_div_ps(rep,ymm2);
                        cv.im = _mm256_div_ps(imp,ymm2);
                        return (cv);
             }
             
#include "GMS_simd_utils.hpp"


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm8c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm256_fmsub_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                             const float * __restrict __ATTR_ALIGN__(32) xim,
                                             const float * __restrict __ATTR_ALIGN__(32) yre,
                                             const float * __restrict __ATTR_ALIGN__(32) yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm256_fmsub_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
              } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm8c4(const __m256 xre,
                                           const __m256 xim,
                                           const __m256 yre,
                                           const __m256 yim,
                                           float * __restrict mre,
                                           float * __restrict min) {

                        register __m256 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre = 0.0f;
                        rep  = _mm256_fmsub_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm8c4(const ymm8c4_t x,
                                           const ymm8c4_t y,
                                           float * __restrict mre,
                                           float * __restrict min) {

                        register __m256 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre = 0.0f;
                        rep  = _mm256_fmsub_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
             }


               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm8c4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        imp  = _mm256_fmadd_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_ps(ymm1,ymm1,
                                               _mm256_mul_ps(ymm3,ymm3));
                        rquot = _mm256_div_ps(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_ps(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                             const float * __restrict __ATTR_ALIGN__(32) xim,
                                             const float * __restrict __ATTR_ALIGN__(32) yre,
                                             const float * __restrict __ATTR_ALIGN__(32) yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        imp  = _mm256_fmadd_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_ps(ymm1,ymm1,
                                               _mm256_mul_ps(ymm3,ymm3));
                        rquot = _mm256_div_ps(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_ps(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm8c4(  const __m256 xre,
                                             const __m256 xim,
                                             const __m256 yre,
                                             const __m256 yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        imp  = _mm256_fmadd_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_ps(yre,yre,
                                               _mm256_mul_ps(yim,yim));
                        rquot = _mm256_div_ps(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_ps(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm8c4(  const ymm8c4_t x,
                                             const ymm8c4_t y,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        imp  = _mm256_fmadd_ps(xim,yre,
                                               _mm256_mul_ps(x.re,y.im));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_ps(y.re,y.re,
                                               _mm256_mul_ps(y.im,y.im));
                        rquot = _mm256_div_ps(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_ps(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm8c4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim ) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp,magc1,magc2,vcmag;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        rep  = _mm256_fmadd_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        _mm256_storeu_ps(&mre[0], _mm256_div_ps(rep,vcmag));
                        _mm256_storeu_ps(&mim[0], _mm256_div_ps(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                              const float * __restrict __ATTR_ALIGN__(32) xim,
                                              const float * __restrict __ATTR_ALIGN__(32) yre,
                                              const float * __restrict __ATTR_ALIGN__(32) yim,
                                              float * __restrict __ATTR_ALIGN__(32) mre,
                                              float * __restrict __ATTR_ALIGN__(32) mim ) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 rep,imp,magc1,magc2,vcmag;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        rep  = _mm256_fmad_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        _mm256_store_ps(&mre[0], _mm256_div_ps(rep,vcmag));
                        _mm256_store_ps(&mim[0], _mm256_div_ps(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm8c4(const __m256 xre,
                                            const __m256 xim,
                                            const __m256 yre,
                                            const __m256 yim,
                                            __m256 * __restrict mre,
                                            __m256 * __restrict mim) {

                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        *mre = _mm256_div_ps(rep,vcmag);
                        *mim = _mm256_div_ps(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm8c4(const ymm8c4_t x,
                                            const ymm8c4_t y,
                                            __m256 * __restrict mre,
                                            __m256 * __restrict mim) {

                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        *mre = _mm256_div_ps(rep,vcmag);
                        *mim = _mm256_div_ps(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cnorm_cprod_ymm8c4(const __m256 xre,
                                                 const __m256 xim,
                                                 const __m256 yre,
                                                 const __m256 yim) {
                                               
                        ymm8c4_t cv;
                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        cv.re = _mm256_div_ps(rep,vcmag);
                        cv.im = _mm256_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cnorm_cprod_ymm8c4(const ymm8c4_t x,
                                                 const ymm8c4_t y) {
                                               
                        ymm8c4_t cv;
                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        cv.re = _mm256_div_ps(rep,vcmag);
                        cv.im = _mm256_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
             


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm8c4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim) {
                      
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        re   = _mm256_fmadd_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                              const float * __restrict __ATTR_ALIGN__(32) xim,
                                              const float * __restrict __ATTR_ALIGN__(32) yre,
                                              const float * __restrict __ATTR_ALIGN__(32) yim,
                                              float * __restrict mre,
                                              float * __restrict mim) {
                      
                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        re   = _mm256_fmadd_ps(ymm0,ymm1,
                                               _mm256_mul_ps(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_ps(ymm2,ymm1,
                                               _mm256_mul_ps(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(rep);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm8c4(const __m256 xre,
                                            const __m256 xim,
                                            const __m256 yre,
                                            const __m256 yim,
                                            float * __restrict mre,
                                            float * __restrict min) {

                        register __m256 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm256_fmadd_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm8c4(const ymm8c4_t x,
                                            const ymm8c4_t y,
                                            float * __restrict mre,
                                            float * __restrict min) {

                        register __m256 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm256_fmadd_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm8c4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        register __m256 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm256_loadu_ps(&xre[0]);
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_loadu_ps(&xim[0]);
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                              const float * __restrict __ATTR_ALIGN__(32) xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        register __m256 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm256_load_ps(&xre[0]);
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_load_ps(&xim[0]);
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm8c4(  const __m256 xre,
                                              const __m256 xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre  = ymm8r4_horizontal_sum(xre);
                        *mre = sre*inv16;
                        sim  = ymm8r4_horizontal_sum(xim);
                        *mim = sim*inv16; 
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm8c4(  const ymm8c4_t x,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre  = ymm8r4_horizontal_sum(x.re);
                        *mre = sre*inv16;
                        sim  = ymm8r4_horizontal_sum(x.im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm8c4_u( const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim ) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 re,im,cvmag;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(ymm0,ymm1,
                                                              _mm256_mul_ps(ymm2,ymm3)));
                        _mm256_storeu_ps(&mre[0], _mm256_div_ps(ymm0,cvmag));
                        _mm256_storeu_ps(&mim[0], _mm256_div_ps(ymm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm8c4_a( const float * __restrict __ATTR_ALIGN__(32) xre,
                                              const float * __restrict __ATTR_ALIGN__(32) xim,
                                              const float * __restrict __ATTR_ALIGN__(32) yre,
                                              const float * __restrict __ATTR_ALIGN__(32) yim,
                                              float * __restrict __ATTR_ALIGN__(32) mre,
                                              float * __restrict __ATTR_ALIGN__(32) mim ) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 re,im,cvmag;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(ymm0,ymm1,
                                                              _mm256_mul_ps(ymm2,ymm3)));
                        _mm256_store_ps(&mre[0], _mm256_div_ps(ymm0,cvmag));
                        _mm256_store_ps(&mim[0], _mm256_div_ps(ymm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm8c4( const __m256 xre,
                                            const __m256 xim,
                                            const __m256 yre,
                                            const __m256 yim,
                                            __m256 * __restrict mre,
                                            __m256 * __restrict mim ) {

                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(xre,yre,
                                                    _mm256_mul_ps(xim,yim)));
                        *mre = _mm256_div_ps(xre,cvmag));
                        *mim =  _mm256_div_ps(xim,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm8c4( const ymm8c4_t x,
                                            const ymm8c4_t y,
                                            __m256 * __restrict mre,
                                            __m256 * __restrict mim ) {

                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(x.re,y.re,
                                                    _mm256_mul_ps(x.im,y.im)));
                        *mre = _mm256_div_ps(x.re,cvmag));
                        *mim =  _mm256_div_ps(x.im,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cnormalize_ymm8c4( const __m256 xre,
                                                 const __m256 xim,
                                                 const __m256 yre,
                                                 const __m256 yim) {
                                            
                        ymm8c4_t cv;
                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(xre,yre,
                                                    _mm256_mul_ps(xim,yim)));
                        cv.re = _mm256_div_ps(xre,cvmag));
                        cv.im =  _mm256_div_ps(xim,cvmag));
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm8c4_t cnormalize_ymm8c4( const ymm8c4_t x,
                                                 const ymm8c4_t y,) {
                                            
                        ymm8c4_t cv;
                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(x.re,y.re,
                                                    _mm256_mul_ps(x.im,y.im)));
                        cv.re = _mm256_div_ps(x.re,cvmag));
                        cv.im =  _mm256_div_ps(x.im,cvmag));
                        return (cv);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm8c4_u( const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 cvmag;
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&yre[0]);
                        ymm2 = _mm256_loadu_ps(&xim[0]);
                        ymm3 = _mm256_loadu_ps(&yim[0]);
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(ymm0,ymm1,
                                                          _mm256_mul_ps(ymm2,ymm3)));
                        _mm256_storeu_ps(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm8c4_a( const float * __restrict __ATTR_ALIGN__(32) xre,
                                              const float * __restrict __ATTR_ALIGN__(32) xim,
                                              const float * __restrict __ATTR_ALIGN__(32) yre,
                                              const float * __restrict __ATTR_ALIGN__(32) yim,
                                              float * __restrict __ATTR_ALIGN__(32) mre) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 cvmag;
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&yre[0]);
                        ymm2 = _mm256_load_ps(&xim[0]);
                        ymm3 = _mm256_load_ps(&yim[0]);
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(ymm0,ymm1,
                                                          _mm256_mul_ps(ymm2,ymm3)));
                        _mm256_store_ps(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm8c4(   const __m256 xre,
                                              const __m256 xim,
                                              const __m256 yre,
                                              const __m256 yim,
                                              __m256 * __restrict  mre) {

                        register __m256 cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(xre,yre,
                                                          _mm256_mul_ps(xim,yim)));
                        *mre = cvmag;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm8c4(   const ymm8c4_t x,
                                              const ymm8c4_t y,
                                              __m256 * __restrict  mre) {

                        register __m256 cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(x.re,y.re,
                                                          _mm256_mul_ps(x.im,y.im)));
                        *mre = cvmag;
             }


#include <complex>


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_2xr4_c4_unroll16x(float * __restrict __ATTR_ALIGN__(32) xre,
                                               float * __restrict __ATTR_ALIGN__(32) xim,
                                               std::complex<float> * __restrict __ATTR_ALIGN__(32) vc,
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


                /*   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_2xr8_c8_unroll16x(double * __restrict __ATTR_ALIGN__(32) xre,
                                               double * __restrict __ATTR_ALIGN__(32) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(32) vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register double re0,im0,re1,im1,re2,im2;
                        register double re3,im3,re4,im4,re5,im5;
                        register double re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (double*)__builtin_assume_aligned(xre,64);
                        xim = (double*)__builtin_assume_aligned(xim,64);
                        vc  = (std::complex<double>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const double re = xre[j];
                               const double im = xim[j];
                               const std::complex<double> tj{re,im};
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

*/
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_c4_2xr4_unroll16x( float * __restrict __ATTR_ALIGN__(32) xre,
                                               float * __restrict __ATTR_ALIGN__(32) xim,
                                               std::complex<float> * __restrict __ATTR_ALIGN__(32) vc,
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


               /*    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_c8_2xr8_unroll16x( double * __restrict __ATTR_ALIGN__(32) xre,
                                                double * __restrict __ATTR_ALIGN__(32) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(32) vc,
                                               const int32_t n) {

                        if(__builtin_expect(0<=n,0)) { return;}
                        register double re0,im0,re1,im1,re2,im2;
                        register double re3,im3,re4,im4,re5,im5;
                        register double re6,im6,re7,im7;
                        int32_t j,m,m1;
                        m = n%16;
#if defined(__INTEL_COMPILER)
                        __assume_aligned(xre,64);
                        __assume_aligned(xim,64);
                        __assume_aligned(c,64);
#elif defined(__GNUC__) && (!defined __INTEL_COMPILER)
                        xre = (double*)__builtin_assume_aligned(xre,64);
                        xim = (double*)__builtin_assume_aligned(xim,64);
                        vc   = (std::complex<double>*)__builtin_assume_aligned(vc,64);
#endif
                        if(n != 0) {
#pragma omp simd aligned(xre:64) aligned(xim:64) aligned(c:64)
                           for(j = 0; j != m; ++j) {
                               const double re = vc[j].real();
                               const double im = vc[j].imag();
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
               */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_ymm8c4_u(const float * __restrict xre,
                                               const float * __restrict xim,
                                               std::complex<float> * __restrict yc) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 ymm4,ymm5;
                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        ymm0 = _mm256_loadu_ps(&xre[0]);
                        ymm1 = _mm256_loadu_ps(&xim[0]);
                        ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                        ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                        ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                        _mm256_storeu_ps(&pyc[0],ymm4);
                        ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                        _mm256_storeu_ps(&pyc[8],ymm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_ymm8c4(  const __m256 xre,
                                               const __m256 xim,
                                               std::complex<float> * __restrict yc) {

                        register __m256 ymm2,ymm3;
                        register __m256 ymm4,ymm5;
                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        ymm2 = _mm256_unpacklo_ps(xre,xim);
                        ymm3 = _mm256_unpackhi_ps(xre,xim);
                        ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                        _mm256_storeu_ps(&pyc[0],ymm4);
                        ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                        _mm256_storeu_ps(&pyc[8],ymm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_ymm8c4_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                               const float * __restrict __ATTR_ALIGN__(32) xim,
                                               std::complex<float> * __restrict __ATTR_ALIGN__(32) yc) {

                        register __m256 ymm0,ymm1,ymm2,ymm3;
                        register __m256 ymm4,ymm5;
                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        ymm0 = _mm256_load_ps(&xre[0]);
                        ymm1 = _mm256_load_ps(&xim[0]);
                        ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                        ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                        ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                        _mm256_store_ps(&pyc[0],ymm4);
                        ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                        _mm256_store_ps(&pyc[8],ymm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void copy_2r4c4_ymm8c4_blocked_a(const float * __restrict __ATTR_ALIGN__(32) xre,
                                                   const float * __restrict __ATTR_ALIGN__(32) xim,
                                                   std::complex<float> * __restrict __ATTR_ALIGN__(32) yc,
                                                   const int32_t n) { // size of array std::complex<float> len elements.

                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        if(n <= 1) {
                            return;
                        }
                        else if(n <= 8) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_load_ps(&xre[0]);
                            ymm1 = _mm256_load_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[8],ymm5);
                        }
                        else if(n <= 16) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_load_ps(&xre[0]);
                            ymm1 = _mm256_load_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[8]);
                            ymm1 = _mm256_load_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[24],ymm5);
                       }
                       else if(n <= 32) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_load_ps(&xre[0]);
                            ymm1 = _mm256_load_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[8]);
                            ymm1 = _mm256_load_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[16]);
                            ymm1 = _mm256_load_ps(&xim[16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[24]);
                            ymm1 = _mm256_load_ps(&xim[24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[56],ymm5);
                       }
                       else if(n <= 64) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_load_ps(&xre[0]);
                            ymm1 = _mm256_load_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[8]);
                            ymm1 = _mm256_load_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[16]);
                            ymm1 = _mm256_load_ps(&xim[16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[24]);
                            ymm1 = _mm256_load_ps(&xim[24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[56],ymm5);
                            ymm0 = _mm256_load_ps(&xre[32]);
                            ymm1 = _mm256_load_ps(&xim[32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[72],ymm5);
                            ymm0 = _mm256_load_ps(&xre[40]);
                            ymm1 = _mm256_load_ps(&xim[40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[88],ymm5);
                            ymm0 = _mm256_load_ps(&xre[48]);
                            ymm1 = _mm256_load_ps(&xim[48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[104],ymm5);
                            ymm0 = _mm256_load_ps(&xre[56]);
                            ymm1 = _mm256_load_ps(&xim[56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[120],ymm5);
                       }
                       else if(n <= 128) {
#if (GMS_MAN_PREFETCH) == 1
                               _mm_prefetch((const char*)&xre[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[16],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[16],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[48],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[48],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[80],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[80],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[112],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[112],_MM_HINT_T0);
#endif                           
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_load_ps(&xre[0]);
                            ymm1 = _mm256_load_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[8]);
                            ymm1 = _mm256_load_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[16]);
                            ymm1 = _mm256_load_ps(&xim[16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[24]);
                            ymm1 = _mm256_load_ps(&xim[24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[56],ymm5);
                            ymm0 = _mm256_load_ps(&xre[32]);
                            ymm1 = _mm256_load_ps(&xim[32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[72],ymm5);
                            ymm0 = _mm256_load_ps(&xre[40]);
                            ymm1 = _mm256_load_ps(&xim[40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[88],ymm5);
                            ymm0 = _mm256_load_ps(&xre[48]);
                            ymm1 = _mm256_load_ps(&xim[48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[104],ymm5);
                            ymm0 = _mm256_load_ps(&xre[56]);
                            ymm1 = _mm256_load_ps(&xim[56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[120],ymm5);
                            ymm0 = _mm256_load_ps(&xre[64]);
                            ymm1 = _mm256_load_ps(&xim[64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[136],ymm5);
                            ymm0 = _mm256_load_ps(&xre[72]);
                            ymm1 = _mm256_load_ps(&xim[72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[152],ymm5);
                            ymm0 = _mm256_load_ps(&xre[80]);
                            ymm1 = _mm256_load_ps(&xim[80]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[168],ymm5);
                            ymm0 = _mm256_load_ps(&xre[88]);
                            ymm1 = _mm256_load_ps(&xim[88]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[184],ymm5);
                            ymm0 = _mm256_load_ps(&xre[96]);
                            ymm1 = _mm256_load_ps(&xim[96]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[200],ymm5);
                            ymm0 = _mm256_load_ps(&xre[104]);
                            ymm1 = _mm256_load_ps(&xim[104]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[216],ymm5);
                            ymm0 = _mm256_load_ps(&xre[112]);
                            ymm1 = _mm256_load_ps(&xim[112]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[232],ymm5);
                            ymm0 = _mm256_load_ps(&xre[120]);
                            ymm1 = _mm256_load_ps(&xim[120]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[248],ymm5);
                       }
                      else if(n > 128) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            int32_t i;
                       for(i = 0; (i+127) < n; i += 128) {
#if (GMS_MAN_PREFETCH) == 1
                            _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+112],_MM_HINT_T0);
#endif                         
                            ymm0 = _mm256_load_ps(&xre[i+0]);
                            ymm1 = _mm256_load_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+8]);
                            ymm1 = _mm256_load_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+16]);
                            ymm1 = _mm256_load_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+24]);
                            ymm1 = _mm256_load_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+56],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+32]);
                            ymm1 = _mm256_load_ps(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+72],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+40]);
                            ymm1 = _mm256_load_ps(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+88],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+48]);
                            ymm1 = _mm256_load_ps(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+104],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+56]);
                            ymm1 = _mm256_load_ps(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+120],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+64]);
                            ymm1 = _mm256_load_ps(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+136],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+72]);
                            ymm1 = _mm256_load_ps(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+152],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+80]);
                            ymm1 = _mm256_load_ps(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+168],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+88]);
                            ymm1 = _mm256_load_ps(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+184],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+96]);
                            ymm1 = _mm256_load_ps(&xim[i+96]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+200],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+104]);
                            ymm1 = _mm256_load_ps(&xim[i+104]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+216],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+112]);
                            ymm1 = _mm256_load_ps(&xim[i+112]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+232],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+120]);
                            ymm1 = _mm256_load_ps(&xim[i+120]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+248],ymm5);
                     }

                       for(; (i+95) < n; i += 96) {
                            ymm0 = _mm256_load_ps(&xre[i+0]);
                            ymm1 = _mm256_load_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+8]);
                            ymm1 = _mm256_load_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+16]);
                            ymm1 = _mm256_load_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+24]);
                            ymm1 = _mm256_load_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+56],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+32]);
                            ymm1 = _mm256_load_ps(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+72],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+40]);
                            ymm1 = _mm256_load_ps(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+88],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+48]);
                            ymm1 = _mm256_load_ps(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+104],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+56]);
                            ymm1 = _mm256_load_ps(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+120],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+64]);
                            ymm1 = _mm256_load_ps(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+136],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+72]);
                            ymm1 = _mm256_load_ps(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+152],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+80]);
                            ymm1 = _mm256_load_ps(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+168],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+88]);
                            ymm1 = _mm256_load_ps(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+184],ymm5);
                     }    

                       for(; (i+79) < n; i += 80) {
                            ymm0 = _mm256_load_ps(&xre[i+0]);
                            ymm1 = _mm256_load_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+8]);
                            ymm1 = _mm256_load_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+16]);
                            ymm1 = _mm256_load_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+24]);
                            ymm1 = _mm256_load_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+56],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+32]);
                            ymm1 = _mm256_load_ps(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+72],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+40]);
                            ymm1 = _mm256_load_ps(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+88],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+48]);
                            ymm1 = _mm256_load_ps(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+104],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+56]);
                            ymm1 = _mm256_load_ps(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+120],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+64]);
                            ymm1 = _mm256_load_ps(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+136],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+72]);
                            ymm1 = _mm256_load_ps(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+152],ymm5);
                     }

                       for(; (i+31) < n; i += 32) {
                            ymm0 = _mm256_load_ps(&xre[i+0]);
                            ymm1 = _mm256_load_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+8]);
                            ymm1 = _mm256_load_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+16]);
                            ymm1 = _mm256_load_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_ps(&xre[i+24]);
                            ymm1 = _mm256_load_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+56],ymm5);
                     }

                     for(; (i+15) < n; i += 16) {
                            ymm0 = _mm256_load_ps(&xre[i+0]);
                            ymm1 = _mm256_load_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+8],ymm5);
                            
                     }

                       for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             pyc[i]         = re;
                             pyc[i+1]       = im;
                      
                       }

                }


 
            }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void copy_2r4c4_ymm8c4_blocked_u(const float * __restrict  xre,
                                                   const float * __restrict  xim,
                                                   std::complex<float> *  yc,
                                                   const int32_t n) { // size of array std::complex<float> len elements.

                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        if(n <= 1) {
                            return;
                        }
                        else if(n <= 8) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_loadu_ps(&xre[0]);
                            ymm1 = _mm256_loadu_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[8],ymm5);
                        }
                        else if(n <= 16) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_loadu_ps(&xre[0]);
                            ymm1 = _mm256_loadu_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[8]);
                            ymm1 = _mm256_loadu_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[24],ymm5);
                       }
                       else if(n <= 32) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_loadu_ps(&xre[0]);
                            ymm1 = _mm256_loadu_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[8]);
                            ymm1 = _mm256_loadu_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[16]);
                            ymm1 = _mm256_loadu_ps(&xim[16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[24]);
                            ymm1 = _mm256_loadu_ps(&xim[24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[56],ymm5);
                       }
                       else if(n <= 64) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_loadu_ps(&xre[0]);
                            ymm1 = _mm256_loadu_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[8]);
                            ymm1 = _mm256_loadu_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[16]);
                            ymm1 = _mm256_loadu_ps(&xim[16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[24]);
                            ymm1 = _mm256_loadu_ps(&xim[24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[56],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[32]);
                            ymm1 = _mm256_loadu_ps(&xim[32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[72],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[40]);
                            ymm1 = _mm256_loadu_ps(&xim[40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[88],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[48]);
                            ymm1 = _mm256_loadu_ps(&xim[48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[104],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[56]);
                            ymm1 = _mm256_loadu_ps(&xim[56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[120],ymm5);
                       }
                       else if(n <= 128) {
#if (GMS_MAN_PREFETCH) == 1
                               _mm_prefetch((const char*)&xre[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[16],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[16],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[48],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[48],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[80],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[80],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[112],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[112],_MM_HINT_T0);
#endif                           
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            ymm0 = _mm256_loadu_ps(&xre[0]);
                            ymm1 = _mm256_loadu_ps(&xim[0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[8]);
                            ymm1 = _mm256_loadu_ps(&xim[8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[16]);
                            ymm1 = _mm256_loadu_ps(&xim[16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[24]);
                            ymm1 = _mm256_loadu_ps(&xim[24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[56],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[32]);
                            ymm1 = _mm256_loadu_ps(&xim[32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[72],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[40]);
                            ymm1 = _mm256_loadu_ps(&xim[40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[88],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[48]);
                            ymm1 = _mm256_loadu_ps(&xim[48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[104],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[56]);
                            ymm1 = _mm256_loadu_ps(&xim[56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[120],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[64]);
                            ymm1 = _mm256_loadu_ps(&xim[64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[136],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[72]);
                            ymm1 = _mm256_loadu_ps(&xim[72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[152],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[80]);
                            ymm1 = _mm256_loadu_ps(&xim[80]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[168],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[88]);
                            ymm1 = _mm256_loadu_ps(&xim[88]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[184],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[96]);
                            ymm1 = _mm256_loadu_ps(&xim[96]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[200],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[104]);
                            ymm1 = _mm256_loadu_ps(&xim[104]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[216],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[112]);
                            ymm1 = _mm256_loadu_ps(&xim[112]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[232],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[120]);
                            ymm1 = _mm256_loadu_ps(&xim[120]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[248],ymm5);
                       }
                      else if(n > 128) {
                            register __m256 ymm0,ymm1,ymm2,ymm3;
                            register __m256 ymm4,ymm5;
                            int32_t i;
                       for(i = 0; (i+127) < n; i += 128) {
#if (GMS_MAN_PREFETCH) == 1
                            _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+16],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+48],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+80],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+112],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+112],_MM_HINT_T0);
#endif                         
                            ymm0 = _mm256_loadu_ps(&xre[i+0]);
                            ymm1 = _mm256_loadu_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+8]);
                            ymm1 = _mm256_loadu_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+16]);
                            ymm1 = _mm256_loadu_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+24]);
                            ymm1 = _mm256_loadu_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+56],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+32]);
                            ymm1 = _mm256_loadu_ps(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+72],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+40]);
                            ymm1 = _mm256_loadu_ps(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+88],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+48]);
                            ymm1 = _mm256_loadu_ps(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+104],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+56]);
                            ymm1 = _mm256_loadu_ps(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+120],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+64]);
                            ymm1 = _mm256_loadu_ps(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+136],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+72]);
                            ymm1 = _mm256_loadu_ps(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+152],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+80]);
                            ymm1 = _mm256_loadu_ps(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+168],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+88]);
                            ymm1 = _mm256_loadu_ps(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+184],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+96]);
                            ymm1 = _mm256_loadu_ps(&xim[i+96]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+200],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+104]);
                            ymm1 = _mm256_loadu_ps(&xim[i+104]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+216],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+112]);
                            ymm1 = _mm256_loadu_ps(&xim[i+112]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+232],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+120]);
                            ymm1 = _mm256_loadu_ps(&xim[i+120]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+248],ymm5);
                     }

                       for(; (i+95) < n; i += 96) {
                            ymm0 = _mm256_loadu_ps(&xre[i+0]);
                            ymm1 = _mm256_loadu_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+8]);
                            ymm1 = _mm256_loadu_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+16]);
                            ymm1 = _mm256_loadu_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+24]);
                            ymm1 = _mm256_loadu_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+56],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+32]);
                            ymm1 = _mm256_loadu_ps(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+72],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+40]);
                            ymm1 = _mm256_loadu_ps(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+88],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+48]);
                            ymm1 = _mm256_loadu_ps(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+104],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+56]);
                            ymm1 = _mm256_loadu_ps(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+120],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+64]);
                            ymm1 = _mm256_loadu_ps(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_ps(&pyc[i+136],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+72]);
                            ymm1 = _mm256_loadu_ps(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+152],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+80]);
                            ymm1 = _mm256_loadu_ps(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+168],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+88]);
                            ymm1 = _mm256_loadu_ps(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+184],ymm5);
                     }    

                       for(; (i+79) < n; i += 80) {
                            ymm0 = _mm256_loadu_ps(&xre[i+0]);
                            ymm1 = _mm256_loadu_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+8]);
                            ymm1 = _mm256_loadu_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+16]);
                            ymm1 = _mm256_loadu_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+24]);
                            ymm1 = _mm256_loadu_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+56],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+32]);
                            ymm1 = _mm256_loadu_ps(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+72],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+40]);
                            ymm1 = _mm256_loadu_ps(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+88],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+48]);
                            ymm1 = _mm256_loadu_ps(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+104],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+56]);
                            ymm1 = _mm256_loadu_ps(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+120],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+64]);
                            ymm1 = _mm256_loadu_ps(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+136],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+72]);
                            ymm1 = _mm256_loadu_ps(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+152],ymm5);
                     }

                       for(; (i+31) < n; i += 32) {
                            ymm0 = _mm256_loadu_ps(&xre[i+0]);
                            ymm1 = _mm256_loadu_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+8]);
                            ymm1 = _mm256_loadu_ps(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+16]);
                            ymm1 = _mm256_loadu_ps(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_ps(&xre[i+24]);
                            ymm1 = _mm256_loadu_ps(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_ps(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+56],ymm5);
                     }

                     for(; (i+15) < n; i += 16) {
                            ymm0 = _mm256_loadu_ps(&xre[i+0]);
                            ymm1 = _mm256_loadu_ps(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_ps(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_ps(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_ps(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_ps(&pyc[i+8],ymm5);
                           
                     }

                       for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             pyc[i]         = re;
                             pyc[i+1]       = im;
                      
                       }

                }


 
            }
   



      } // math


} // gms















#endif /*__GMS_COMPLEX_ZMM16R4_HPP__*/
