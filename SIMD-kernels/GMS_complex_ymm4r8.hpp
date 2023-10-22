
#ifndef __GMS_COMPLEX_YMM4R8_HPP__
#define __GMS_COMPLEX_YMM4R8_HPP__ 181020231557


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

    const unsigned int GMS_COMPLEX_YMM4R8_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_YMM4R8_MINOR = 0U;
    const unsigned int GMS_COMPLEX_YMM4R8_MICRO = 0U;
    const unsigned int GMS_COMPLEX_YMM4R8_FULLVER =
      1000U*GMS_COMPLEX_YMM4R8_MAJOR+
      100U*GMS_COMPLEX_YMM4R8_MINOR+
      10U*GMS_COMPLEX_YMM4R8_MICRO;
    const char * const GMS_COMPLEX_YMM4R8_CREATION_DATE = "20-10-2023 09:13 AM +00200 (FRI 20 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_YMM4R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_YMM4R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_YMM4R8_DESCRIPTION   = "AVX/AVX2 optimized complex number implementation."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.hpp"

namespace  gms {


       namespace math {
       
       
                   struct __ATTR_ALIGN__(32) ymm4c8_t {
                   
                          __m256d re;
                          __m256d im;
                   };
                   
                   

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_pd(&xre[0]);
                        ymm1  = _mm256_loadu_pd(&yre[0]);
                        _mm256_storeu_pd(&zre[0], _mm256_add_pd(ymm0,ymm1));
                        ymm2  = _mm256_loadu_pd(&xim[0]);
                        ymm3  = _mm256_loadu_pd(&yim[0]);
                        _mm256_storeu_pd(&zim[0], _mm256_add_pd(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       const double * __restrict __ATTR_ALIGN__(32) yre,
                                       const double * __restrict __ATTR_ALIGN__(32) yim,
                                       double *       __restrict __ATTR_ALIGN__(32) zre,
                                       double *       __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_pd(&xre[0]);
                        ymm1  = _mm256_load_pd(&yre[0]);
                        _mm256_store_pd(&zre[0], _mm256_add_pd(ymm0,ymm1));
                        ymm2  = _mm256_load_pd(&xim[0]);
                        ymm3  = _mm256_load_pd(&yim[0]);
                        _mm256_store_pd(&zim[0], _mm256_add_pd(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d yre,
                                     const __m256d yim,
                                     __m256d * __restrict zre,
                                     __m256d * __restrict zim) {
                     
                        register __m256d ymm0,ymm1;
                        ymm0  = _mm256_add_pd(xre,yre);
                        *zre  = ymm0;
                        ymm1  = _mm256_add_pd(xim,yim);
                        *zim  = ymm1;
                }
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cadd_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d yre,
                                          const __m256d yim) {
                                     
                        ymm4c8_t cv;
                        register __m256d ymm0,ymm1;
                        ymm0   = _mm256_add_pd(xre,yre);
                        cv.re  = ymm0;
                        ymm1   = _mm256_add_pd(xim,yim);
                        cv.im  = ymm1;  
                        return (cv);            
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cadd_ymm4c8(const ymm4c8_t x,
                                          const ymm4c8_t y) {
                                     
                        ymm4c8_t cv;
                        register __m256d ymm0,ymm1;
                        ymm0   = _mm256_add_pd(x.re,y.re);
                        cv.re  = ymm0;
                        ymm1   = _mm256_add_pd(x.im,y.im);
                        cv.im  = ymm1;  
                        return (cv);            
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d s,
                                     __m256d * __restrict     zre,
                                     __m256d * __restrict     zim) {

                        *zre = _mm256_add_pd(xre,s);
                        *zim = _mm256_setzero_pd();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ymm4c8_t cadd_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d s) {
                      
                      ymm4c8_t cv;
                      cv.re =  _mm256_add_pd(xre,s);
                      cv.im =  _mm256_setzero_pd();
                      return (cv);                       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ymm4c8_t cadd_ymm4c8(const ymm4c8_t x,
                                          const __m256d s) {
                      
                      ymm4c8_t cv;
                      cv.re =  _mm256_add_pd(x.re,s);
                      cv.im =  _mm256_setzero_pd();
                      return (cv);                       
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm4c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_pd(&xre[0]);
                        ymm1  = _mm256_loadu_pd(&xim[0]);
                        ymm2  = _mm256_loadu_pd(&zre[0]);
                        ymm3  = _mm256_loadu_pd(&zim[0])
                        _mm256_storeu_pd(&zre[0], _mm256_add_pd(ymm2,ymm0));
                        _mm256_storeu_pd(&zim[0], _mm256_add_pd(ymm3,ymm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_ymm4c8_aip(const double * __restrict __ATTR_ALIGN__(32) xre,
                                         const double * __restrict __ATTR_ALIGN__(32) xim,
                                         double *       __restrict __ATTR_ALIGN__(32) zre,
                                         double *       __restrict __ATTR_ALIGN__(32) zim) {
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_pd(&xre[0]);
                        ymm1  = _mm256_load_pd(&xim[0]);
                        ymm2  = _mm256_load_pd(&zre[0]);
                        ymm3  = _mm256_load_pd(&zim[0])
                        _mm256_store_pd(&zre[0], _mm256_add_pd(ymm2,ymm0));
                        _mm256_store_pd(&zim[0], _mm256_add_pd(ymm3,ymm1));
              }


                ////////////////////////////////////////////////////////////////////

                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_pd(&xre[0]);
                        ymm1  = _mm256_loadu_pd(&yre[0]);
                        _mm256_storeu_pd(&zre[0], _mm256_sub_pd(ymm0,ymm1));
                        ymm2  = _mm256_loadu_pd(&xim[0]);
                        ymm3  = _mm256_loadu_pd(&yim[0]);
                        _mm256_storeu_pd(&zim[0], _mm256_sub_pd(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm4c8_a( const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       const double * __restrict __ATTR_ALIGN__(32) yre,
                                       const double * __restrict __ATTR_ALIGN__(32) yim,
                                       double *       __restrict __ATTR_ALIGN__(32) zre,
                                       double *       __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_pd(&xre[0]);
                        ymm1  = _mm256_load_pd(&yre[0]);
                        _mm256_store_pd(&zre[0], _mm256_sub_pd(ymm0,ymm1));
                        ymm2  = _mm256_load_pd(&xim[0]);
                        ymm3  = _mm256_load_pd(&yim[0]);
                        _mm256_store_pd(&zim[0], _mm256_sub_pd(ymm2,ymm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d yre,
                                     const __m256d yim,
                                     __m256d * __restrict     zre,
                                     __m256d * __restrict     zim) {
                     
                        register __m256d ymm0,ymm1;
                        ymm0  = _mm256_sub_pd(xre,yre);
                        *zre  = ymm0;
                        ymm1  = _mm256_sub_pd(xim,yim);
                        *zim  = ymm1;
                }
                
                
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t csub_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d yre,
                                          const __m256d yim) {
                                    
                        ymm4c8_t cv;
                        register __m256d ymm0,ymm1;
                        ymm0  = _mm256_sub_pd(xre,yre);
                        cv.re  = ymm0;
                        ymm1  = _mm256_sub_pd(xim,yim);
                        cv.im  = ymm1;
                        return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t csub_ymm4c8(const ymm4c8_t x,
                                          const ymm4c8_t y) {
                                    
                        ymm4c8_t cv;
                        register __m256d ymm0,ymm1;
                        ymm0  = _mm256_sub_pd(x.re,y.re);
                        cv.re  = ymm0;
                        ymm1  = _mm256_sub_pd(x.im,y.im);
                        cv.im  = ymm1;
                        return (cv);
                }
                
                


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d s,
                                     __m256d * __restrict     zre,
                                     __m256d * __restrict     zim) {

                        *zre = _mm256_sub_pd(xre,s);
                        *zim = _mm256_setzero_pd();
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t csub_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d s) {
                                    
                        ymm4c8_t cv;
                        cv.re = _mm256_sub_pd(xre,s);
                        cv.im = _mm256_setzero_pd();
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t csub_ymm4c8(const ymm4c8_t x,
                                          const __m256d s) {
                                    
                        ymm4c8_t cv;
                        cv.re = _mm256_sub_pd(x.re,s);
                        cv.im = _mm256_setzero_pd();
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm4c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_pd(&xre[0]);
                        ymm1  = _mm256_loadu_pd(&xim[0]);
                        ymm2  = _mm256_loadu_pd(&zre[0]);
                        ymm3  = _mm256_loadu_pd(&zim[0])
                        _mm256_storeu_pd(&zre[0], _mm256_sub_pd(ymm2,ymm0));
                        _mm256_storeu_pd(&zim[0], _mm256_sub_pd(ymm3,ymm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_ymm4c8_aip(const double * __restrict __ATTR_ALIGN__(32) xre,
                                         const double * __restrict __ATTR_ALIGN__(32) xim,
                                         double *       __restrict __ATTR_ALIGN__(32) zre,
                                         double *       __restrict __ATTR_ALIGN__(32) zim) {
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_pd(&xre[0]);
                        ymm1  = _mm256_load_pd(&xim[0]);
                        ymm2  = _mm256_load_pd(&zre[0]);
                        ymm3  = _mm256_load_pd(&zim[0])
                        _mm256_store_pd(&zre[0], _mm256_sub_pd(ymm2,ymm0));
                        _mm256_store_pd(&zim[0], _mm256_sub_pd(ymm3,ymm1));
              }


               ////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                           register __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_loadu_pd(&xre[0]);
                           ymm1  = _mm256_loadu_pd(&yre[0]);
                           ymm2  = _mm256_loadu_pd(&xim[0]);
                           ymm3  = _mm256_loadu_pd(&yim[0]);
                           ymm4  = _mm256_sub_pd(_mm256_mul_pd(ymm0,ymm1),
                                                                        _mm256_mul_pd(ymm2,ymm3));
                           _mm256_storeu_pd(&zre[0], ymm4);
                           ymm5  = _mm256_mul_pd(_mm256_mul_pd(ymm2,ymm1),
                                                                        _mm256_mul_pd(ymm0,ymm3));
                           _mm256_storeu_pd(&zim[0], ymm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       const double * __restrict __ATTR_ALIGN__(32) yre,
                                       const double * __restrict __ATTR_ALIGN__(32) yim,
                                       double *       __restrict __ATTR_ALIGN__(32) zre,
                                       double *       __restrict __ATTR_ALIGN__(32) zim) {

                           register __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_load_pd(&xre[0]);
                           ymm1  = _mm256_load_pd(&yre[0]);
                           ymm2  = _mm256_load_pd(&xim[0]);
                           ymm3  = _mm256_load_pd(&yim[0]);
                           ymm4  = _mm256_sub_pd(_mm256_mul_pd(ymm0,ymm1),
                                                                        _mm256_mul_pd(ymm2,ymm3));
                           _mm256_store_pd(&zre[0], ymm4);
                           ymm5  = _mm256_mul_pd(_mm256_mul_pd(ymm2,ymm1),
                                                                        _mm256_mul_pd(ymm0,ymm3));
                           _mm256_store_pd(&zim[0], ymm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d yre,
                                     const __m256d yim,
                                     __m256d * __restrict     zre,
                                     __m256d * __restrict     zim) {

                         register __m256d ymm0,ymm1;
                         ymm0 = _mm256_sub_pd(_mm256_mul_pd(xre,yre),
                                              _mm256_mul_pd(xim,yim));
                         *zre  = ymm0;
                         ymm1 = _mm256_mul_pd(_mm256_mul_pd(xim,yre),
                                              _mm256_mul_pd(xre,yim));
                         *zim  = ymm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cmul_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d yre,
                                          const __m256d yim) {
                                     
                         ymm4c8_t cv
                         register __m256d ymm0,ymm1;
                         ymm0 = _mm256_sub_pd(_mm256_mul_pd(xre,yre),
                                              _mm256_mul_pd(xim,yim));
                         cv.re  = ymm0;
                         ymm1 = _mm256_mul_pd(_mm256_mul_pd(xim,yre),
                                              _mm256_mul_pd(xre,yim));
                         cv.im  = ymm1;
                         return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cmul_ymm4c8(const ymm4c8_t x,
                                          const ymm4c8_t y) {
                                     
                         ymm4c8_t cv
                         register __m256d ymm0,ymm1;
                         ymm0 = _mm256_sub_pd(_mm256_mul_pd(x.re,y.re),
                                              _mm256_mul_pd(x.im,y.im));
                         cv.re  = ymm0;
                         ymm1 = _mm256_mul_pd(_mm256_mul_pd(x.im,y.re),
                                              _mm256_mul_pd(x.re,y.im));
                         cv.im  = ymm1;
                         return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d s,
                                     __m256d * __restrict   zre,
                                     __m256d * __restrict   zim) {

                        *zre = _mm256_mul_pd(xre,s);
                        *zim = _mm256_mul_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cmul_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d s) {
                                     
                        ymm4c8_t cv;
                        cv.re = _mm256_mul_pd(xre,s);
                        cv.im = _mm256_mul_pd(xim,s);
                        return (cv);
               }
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cmul_ymm4c8(const ymm4c8_t x,
                                          const __m256d s) {
                                     
                        ymm4c8_t cv;
                        cv.re = _mm256_mul_pd(x.re,s);
                        cv.im = _mm256_mul_pd(x.im,s);
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm4c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {

                           register __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_loadu_pd(&xre[0]);
                           ymm1  = _mm256_loadu_pd(&zre[0]);
                           ymm2  = _mm256_loadu_pd(&xim[0]);
                           ymm3  = _mm256_loadu_pd(&zim[0]);
                           ymm4  = _mm256_sub_pd(_mm256_mul_pd(ymm0,ymm1),
                                                 _mm256_mul_pd(ymm2,ymm3));
                           _mm256_storeu_pd(&zre[0], ymm4);
                           ymm5  = _mm256_mul_pd(_mm256_mul_pd(ymm2,ymm1),
                                                 _mm256_mul_pd(ymm0,ymm3));
                           _mm256_storeu_pd(&zim[0], ymm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_ymm4c8_aip(const double * __restrict __ATTR_ALIGN__(32) xre,
                                         const double * __restrict __ATTR_ALIGN__(32) xim,
                                         double *       __restrict __ATTR_ALIGN__(32) zre,
                                         double *       __restrict __ATTR_ALIGN__(32) zim) {

                           register __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5;
                           ymm0  = _mm256_load_pd(&xre[0]);
                           ymm1  = _mm256_load_pd(&zre[0]);
                           ymm2  = _mm256_load_pd(&xim[0]);
                           ymm3  = _mm256_load_pd(&zim[0]);
                           ymm4  = _mm256_sub_pd(_mm256_mul_pd(ymm0,ymm1),
                                                 _mm256_mul_pd(ymm2,ymm3));
                           _mm256_store_pd(&zre[0], ymm4);
                           ymm5  = _mm256_mul_pd(_mm256_mul_pd(ymm2,ymm1),
                                                 _mm256_mul_pd(ymm0,ymm3));
                           _mm256_store_pd(&zim[0], ymm5);
               }

                 ////////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3; 
                        register __m256d ymm4,ymm5,ymm6;
                        ymm0  = _mm256_loadu_pd(&xre[0]); //a
                        ymm1  = _mm256_loadu_pd(&yim[0]); //d
                        ymm2  = _mm256_loadu_pd(&xim[0]); //b
                        ymm3  = _mm256_loadu_pd(&yre[0]); //c
                        ymm4  = _mm256_fmadd_pd(ymm0,ymm3,
                                                _mm256_mul_pd(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_pd(ymm2,ymm3,
                                                _mm256_mul_pd(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_pd(ymm3,ymm3),
                                                _mm256_mul_pd(ymm1,ymm1));
                        _mm256_storeu_pd(&zre[0], _mm256_div_pd(ymm4,ymm6));
                        _mm256_storeu_pd(&zim[0], _mm256_div_pd(ymm5,ymm6));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8_a(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3; 
                        register __m256d ymm4,ymm5,ymm6;
                        ymm0  = _mm256_load_pd(&xre[0]); //a
                        ymm1  = _mm256_load_pd(&yim[0]); //d
                        ymm2  = _mm256_load_pd(&xim[0]); //b
                        ymm3  = _mm256_load_pd(&yre[0]); //c
                        ymm4  = _mm256_fmadd_pd(ymm0,ymm3,
                                                _mm256_mul_pd(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_pd(ymm2,ymm3,
                                                _mm256_mul_pd(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_pd(ymm3,ymm3,
                                                _mm256_mul_pd(ymm1,ymm1));
                        _mm256_store_pd(&zre[0], _mm256_div_pd(ymm4,ymm6));
                        _mm256_store_pd(&zim[0], _mm256_div_pd(ymm5,ymm6));
                }
              
                                         
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d yre,
                                     const __m256d yim,
                                     __m256d * __restrict zre,
                                     __m256d * __restrict zim) {

                      register __m256d ymm0,ymm1,ymm2;
                      ymm0 = _mm256_fmadd_pd(xre,yre,
                                           _mm256_mul_pd(xim,yim));
                      ymm1 = _mm256_fmsub_pd(xim,yre,
                                           _mm256_mul_pd(xre,yim));
                      ymm2 = _mm256_fmadd_pd(ymm3,ymm3,
                                           _mm256_mul_pd(ymm1,ymm1));
                      *zre  = _mm256_div_pd(ymm0,ymm2);
                      *zim  = _mm256_div_pd(ymm1,ymm2);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d yre,
                                          const __m256d yim) {
                                     
                      ymm4c8_t
                      register __m256d ymm0,ymm1,ymm2;
                      ymm0 = _mm256_fmadd_pd(xre,yre,
                                           _mm256_mul_pd(xim,yim));
                      ymm1 = _mm256_fmsub_pd(x.im,y.re,
                                           _mm256_mul_pd(xre,yim));
                      ymm2 = _mm256_fmadd_pd(ymm3,ymm3,
                                           _mm256_mul_pd(ymm1,ymm1));
                      cv.re  = _mm256_div_pd(ymm0,ymm2);
                      cv.im  = _mm256_div_pd(ymm1,ymm2);
                      return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_ymm4c8(const ymm4c8_t x,
                                          const ymm4c8_t y) {
                                     
                      ymm4c8_t
                      register __m256d ymm0,ymm1,ymm2;
                      ymm0 = _mm256_fmadd_pd(x.re,y.re,
                                           _mm256_mul_pd(x.im,y.im));
                      ymm1 = _mm256_fmsub_pd(x.im,y.re,
                                           _mm256_mul_pd(x.re,y.im));
                      ymm2 = _mm256_fmadd_pd(ymm3,ymm3,
                                           _mm256_mul_pd(ymm1,ymm1));
                      cv.re  = _mm256_div_pd(ymm0,ymm2);
                      cv.im  = _mm256_div_pd(ymm1,ymm2);
                      return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     const __m256d s,
                                     __m256d * __restrict zre,
                                     __m256d * __restrict zim) {

                        *zre = _mm256_div_pd(xre,s);
                        *zim = _mm256_div_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_ymm4c8(const __m256d xre,
                                          const __m256d xim,
                                          const __m256d s) {
                                     
                         ymm4c8_t cv;
                         cv.re = _mm256_div_pd(xre,s);
                         cv.im = _mm256_div_pd(xim,s);
                         return (cv);
               }
               
               
                 
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_ymm4c8(const ymm4c8_t x,
                                          const __m256d s) {
                                     
                         ymm4c8_t cv;
                         cv.re = _mm256_div_pd(x.re,s);
                         cv.im = _mm256_div_pd(x.im,s);
                         return (cv);
               }
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8_s(const __m256d s,
                                       const __m256d xre,
                                       const __m256d xim,
                                       __m256d * __restrict zre,
                                       __m256d * __restrict zim) {
                        
                        register __m256d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm256_setzero_pd();
                        cdiv_ymm4c8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_ymm4c8_s(const __m256d s,
                                            const __m256d xre,
                                            const __m256d xim) {
                                       
                        ymm4c8_t cv;
                        register __m256d t0r,t0i;
                        t0r = s;
                        t0i = _mm256_setzero_pd();
                        cdiv_ymm4c8(t0r,t0i,xre,xim,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_ymm4c8_s(const __m256d s,
                                            const ymm4c8_t x) {
                                       
                        ymm4c8_t cv;
                        register __m256d t0r,t0i;
                        t0r = s;
                        t0i = _mm256_setzero_pd();
                        cdiv_ymm4c8(t0r,t0i,x.re,x.im,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm4c8_s(const __m256d s,
                                             const __m256d xre,
                                             const __m256d xim,
                                             __m256d * __restrict zre,
                                             __m256d * __restrict zim) {
                                             
                        register __m256d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm256_setzero_pd(); 
                        cdiv_smith_ymm4c8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                   
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_smith_ymm4c8_s(const __m256d s,
                                                  const __m256d xre,
                                                  const __m256d xim) {
                                             
                        ymm4c8_t cv;                    
                        register __m256d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm256_setzero_pd(); 
                        cdiv_smith_ymm4c8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_smith_ymm4c8_s(const __m256d s,
                                                  const ymm4c8_t x) {
                                             
                        ymm4c8_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm256_setzero_pd(); 
                        cv = cdiv_smith_ymm4c8(t0,x);
                        return (cv);                 
                 }
                 
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3; 
                        register __m256d ymm4,ymm5,ymm6;
                        ymm0  = _mm256_loadu_pd(&xre[0]); //a
                        ymm1  = _mm256_loadu_pd(&zim[0]); //d
                        ymm2  = _mm256_loadu_pd(&xim[0]); //b
                        ymm3  = _mm256_loadu_pd(&zre[0]); //c
                        ymm4  = _mm256_fmadd_pd(ymm0,ymm3,
                                                _mm256_mul_pd(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_pd(ymm2,ymm3,
                                                _mm256_mul_pd(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_pd(ymm3,ymm3,
                                                _mm256_mul_pd(ymm1,ymm1));
                        _mm256_storeu_pd(&zre[0], _mm256_div_pd(ymm4,ymm6));
                        _mm256_storeu_pd(&zim[0], _mm256_div_pd(ymm5,ymm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_ymm4c8_aip(const double * __restrict __ATTR_ALIGN__(32) xre,
                                         const double * __restrict __ATTR_ALIGN__(32) xim,
                                         double *       __restrict __ATTR_ALIGN__(32) zre,
                                         double *       __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3; 
                        register __m256d ymm4,ymm5,ymm6;
                        ymm0  = _mm256_load_pd(&xre[0]); //a
                        ymm1  = _mm256_load_pd(&zim[0]); //d
                        ymm2  = _mm256_load_pd(&xim[0]); //b
                        ymm3  = _mm256_load_pd(&zre[0]); //c
                        ymm4  = _mm256_fmadd_pd(ymm0,ymm3,
                                                _mm256_mul_pd(ymm2,ymm1));
                        ymm5  = _mm256_fmsub_pd(ymm2,ymm3,
                                                _mm256_mul_pd(ymm0,ymm1));
                        ymm6  = _mm256_fmadd_pd(ymm3,ymm3,
                                                _mm256_mul_pd(ymm1,ymm1));
                        _mm256_store_pd(&zre[0], _mm256_div_pd(ymm4,ymm6));
                        _mm256_store_pd(&zim[0], _mm256_div_pd(ymm5,ymm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm4c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d r,den;
                        __mmask8 m = 0x0;
                        ymm0 = _mm256_loadu_pd(&yre[0]); // c
                        ymm1 = _mm256_loadu_pd(&yim[0]); // d
                        ymm2 = _mm256_loadu_pd(&xre[0]); // a
                        ymm3 = _mm256_loadu_pd(&xim[0]); // b
                        m    = _mm256_cmp_pd_mask(_mm256_abs_pd(ymm0),
                                                  _mm256_abs_pd(ymm1),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_pd(m,_mm256_div_pd(ymm0,ymm1),
                                                      _mm256_div_pd(ymm1,ymm0)); // r
                        den  = _mm256_mask_blend_pd(m,_mm256_fmadd_pd(r,ymm0,ymm1),
                                                      _mm256_fmadd_pd(r,ymm1,ymm0));
                        _mm256_storeu_pd(&zre[0], _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmadd_pd(ymm2,r,ymm3),den),
                                                _mm256_div_pd(_mm256_fmadd_pd(ymm3,r,ymm2),den)));
                        _mm256_storeu_pd(&zim[0], _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmsub_pd(ymm3,r,ymm2),den),
                                                _mm256_div_pd(_mm256_sub_pd(ymm3,_mm256_mul_pd(ymm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm4c8_a(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d r,den;
                        __mmask8 m = 0x0;
                        ymm0 = _mm256_load_pd(&yre[0]); // c
                        ymm1 = _mm256_load_pd(&yim[0]); // d
                        ymm2 = _mm256_load_pd(&xre[0]); // a
                        ymm3 = _mm256_load_pd(&xim[0]); // b
                        m    = _mm256_cmp_pd_mask(_mm256_abs_pd(ymm0),
                                                  _mm256_abs_pd(ymm1),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_pd(m,_mm256_div_pd(ymm0,ymm1),
                                                      _mm256_div_pd(ymm1,ymm0)); // r
                        den  = _mm256_mask_blend_pd(m,_mm256_fmadd_pd(r,ymm0,ymm1),
                                                      _mm256_fmadd_pd(r,ymm1,ymm0));
                        _mm256_storeu_pd(&zre[0], _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmadd_pd(ymm2,r,ymm3),den),
                                                _mm256_div_pd(_mm256_fmadd_pd(ymm3,r,ymm2),den)));
                        _mm256_storeu_pd(&zim[0], _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmsub_pd(ymm3,r,ymm2),den),
                                                _mm256_div_pd(_mm256_sub_pd(ymm3,_mm256_mul_pd(ymm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_ymm4c8(const __m256d xre,
                                           const __m256d xim,
                                           const __m256d yre,
                                           const __m256d yim,
                                           __m256d * __restrict zre,
                                           __m256d * __restrict zim) {

                        register __m256d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_pd_mask(_mm256_abs_pd(yre),
                                                  _mm256_abs_pd(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_pd(m,_mm256_div_pd(yre,yim),
                                                      _mm256_div_pd(yim,yre)); // r
                        den  = _mm256_mask_blend_pd(m,_mm256_fmadd_pd(r,yre,yim),
                                                      _mm256_fmadd_pd(r,yim,yre));
                        *zre  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmadd_pd(xre,r,xim),den),
                                                _mm256_div_pd(_mm256_fmadd_pd(xim,r,xre),den));
                        *zim  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmsub_pd(xim,r,xre),den),
                                                _mm256_div_pd(_mm256_sub_pd(xim,_mm256_mul_pd(xre,r)),den)));
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_smith_ymm4c8(const __m256d xre,
                                                const __m256d xim,
                                                const __m256d yre,
                                                const __m256d yim) {
                                           
                        ymm4c8_t cv
                        register __m256d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_pd_mask(_mm256_abs_pd(yre),
                                                  _mm256_abs_pd(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_pd(m,_mm256_div_pd(yre,yim),
                                                      _mm256_div_pd(yim,yre)); // r
                        den  = _mm256_mask_blend_pd(m,_mm256_fmadd_pd(r,yre,yim),
                                                      _mm256_fmadd_pd(r,yim,yre));
                        cv.re  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmadd_pd(xre,r,xim),den),
                                                _mm256_div_pd(_mm256_fmadd_pd(xim,r,xre),den));
                        cv.im  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmsub_pd(xim,r,xre),den),
                                                _mm256_div_pd(_mm256_sub_pd(xim,_mm256_mul_pd(xre,r)),den)));
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cdiv_smith_ymm4c8(const ymm4c8_t x,
                                                const ymm4c8_t y) {
                                           
                        ymm4c8_t cv
                        register __m256d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_pd_mask(_mm256_abs_pd(y.re),
                                                  _mm256_abs_pd(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_pd(m,_mm256_div_pd(y.re,y.im),
                                                      _mm256_div_pd(y.im,y.re)); // r
                        den  = _mm256_mask_blend_pd(m,_mm256_fmadd_pd(r,y.re,y.im),
                                                      _mm256_fmadd_pd(r,y.im,y.re));
                        cv.re  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmadd_pd(x.re,r,x.im),den),
                                                _mm256_div_pd(_mm256_fmadd_pd(x.im,r,x.re),den));
                        cv.im  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmsub_pd(x.im,r,x.re),den),
                                                _mm256_div_pd(_mm256_sub_pd(x.im,_mm256_mul_pd(x.re,r)),den)));
                        return (cv);
               }





                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_ymm4c8_u(const double * __restrict re,
                                       const double * __restrict im,
                                       double * __restrict  cabs) {

                        register __m256d ymm0,ymm1,ymm2,ymm3,ymm4;
                        ymm0  = _mm256_loadu_pd(&re[0]);
                        ymm1  = _mm256_mul_pd(ymm0,ymm0);
                        ymm2  = _mm256_loadu_pd(&im[0]);
                        ymm3  = _mm256_mul_pd(ymm2,ymm2);
                        ymm4  = _mm256_sqrt_pd(_mm256_add_pd(ymm1,ymm3));
                        _mm256_storeu_pd(&cabs[0],ymm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) re,
                                       const double * __restrict __ATTR_ALIGN__(32) im,
                                       double * __restrict  __ATTR_ALIGN__(32) cabs) {

                        register __m256d ymm0,ymm1,ymm2,ymm3,ymm4;
                        ymm0  = _mm256_load_pd(&re[0]);
                        ymm1  = _mm256_mul_pd(ymm0,ymm0);
                        ymm2  = _mm256_load_pd(&im[0]);
                        ymm3  = _mm256_mul_pd(ymm2,ymm2);
                        ymm4  = _mm256_sqrt_pd(_mm256_add_pd(ymm1,ymm3));
                        _mm256_store_pd(&cabs[0],ymm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256d cabs_ymm4c8(const __m256d re,
                                       const __m256d im) {

                        register __m256d ymm0,ymm1,cabs;
                        ymm0 = _mm256_mul_pd(re,re);
                        ymm1 = _mm256_mul_pd(im,im);
                        cabs = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        return (cabs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256d cabs_ymm4c8(const ymm4c8_t x) {

                        register __m256d ymm0,ymm1,cabs;
                        ymm0 = _mm256_mul_pd(x.re,x.re);
                        ymm1 = _mm256_mul_pd(x.im,x.im);
                        cabs = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        return (cabs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_ymm4c8_u(const double * __restrict re,
                                       const double * __restrict im,
                                       double * __restrict  carg) {

                        register __m256d ymm0,ymm1;
                        ymm0 = _mm256_loadu_pd(&re[0]);
                        ymm1 = _mm256_loadu_pd(&im[0]);
                        _mm256_storeu_pd(&carg[0],_mm256_atan2_pd(ymm0,ymm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) re,
                                       const double * __restrict __ATTR_ALIGN__(32) im,
                                       double * __restrict  __ATTR_ALIGN__(32) carg) {

                        register __m256d ymm0,ymm1;
                        ymm0 = _mm256_load_pd(&re[0]);
                        ymm1 = _mm256_load_pd(&im[0]);
                        _mm256_store_pd(&carg[0],_mm256_atan2_pd(ymm0,ymm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256d carg_ymm4c8(const __m256d re,
                                       const __m256d im) {

                       register __m256d carg;
                       carg = _mm256_atan2_pd(re,im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m256d carg_ymm4c8(ymm4c8_t x) {

                       register __m256d carg;
                       carg = _mm256_atan2_pd(x.re,x.im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_ymm4c8(const __m256d re,
	                             const __m256d im,
	                             __m256d * __restrict clogr,
	                             __m256d * __restrict clogi) {
	                
	                register __m256d t1,t2,ln;
	                t1  = cabs_ymm4c8(re,im);
	                t2  = carg_ymm4c8(re,im);
	                ln  = _mm256_log_pd(t1);
	                *clogr = ln;
	                *clogi = t2;                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) pre,
	                               const double * __restrict __ATTR_ALIGN__(32) pim,
	                               double * __restrict clogr,
	                               double * __restrict clogi) {
	                
	                register __m256d re = _mm256_load_pd(&pre[0]);
	                register __m256d im = _mm256_load_pd(&pim[0]);
	                register __m256d t1,t2,ln;
	                t1  = cabs_ymm4c8(re,im);
	                t2  = carg_ymm4c8(re,im);
	                ln  = _mm256_log_pd(t1);
	                _mm256_store_pd(&clogr[0] ,ln);
	                _mm256_store_pd(&clogi[0] ,t2);                    
	        }
	        
	        
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void clog_ymm4c8_u(const double * __restrict  pre,
	                               const double * __restrict  pim,
	                               double * __restrict clogr,
	                               double * __restrict clogi) {
	                
	                register __m256d re = _mm256_loadu_pd(&pre[0]);
	                register __m256d im = _mm256_loadu_pd(&pim[0]);
	                register __m256d t1,t2,ln;
	                t1  = cabs_ymm4c8(re,im);
	                t2  = carg_ymm4c8(re,im);
	                ln  = _mm256_log_pd(t1);
	                _mm256_storeu_pd(&clogr[0] ,ln);
	                _mm256_storeu_pd(&clogi[0] ,t2);                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ymm4c8_t clog_ymm4c8(const ymm4c8_t x){
	                                  
	                ymm4c8_t clog;                           
	                register __m256d t1,t2,ln;
	                t1  = cabs_ymm4c8(x.re,x.im);
	                t2  = carg_ymm4c8(x.re,x.im);
	                ln  = _mm256_log_pd(t1);
	                clog.re = ln;
	                clog.im = t2;
	                return (clog);                    
	        }

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm4c8_u(double * __restrict re,
                                        double * __restrict im) {

                        register __m256d c;
                        c = negate_ymm8r4(_mm256_loadu_pd(&im[0]));
                        _mm256_storeu_pd(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm4c8_a(double * __restrict __ATTR_ALIGN__(32) re,
                                        double * __restrict __ATTR_ALIGN__(32) im) {
                                        
                        register __m256d c;
                        c = negate_ymm8r4(_mm256_load_pd(&im[0]));
                        _mm256_store_pd(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm4c8(__m256d * __restrict re,
                                      __m256d * __restrict im) {
                         
                        register __m256d c;              
                        c = negate_ymm8r4(*im);
                        *im = c;
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_ymm4c8_v2(const __m256d xre,
                                         const __m256d xim,
                                         __m256d * __restrict yre,
                                         __m256d * __restrict yim) {
                         
                        //register __m256d c;              
                        //c = negate_ymm4c8(*im);
                        //*im = c;
                        *yre = xre; 
                        *yim = negate_ymm8r4(xim);
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cconj_ymm4c8_v2(const __m256d xre,
                                              const __m256d xim) {                                              
                         
                        //register __m256d c;              
                        //c = negate_ymm4c8(*im);
                        //*im = c;
                        ymm4c8_t cv;
                        cv.re = xre; 
                        cv.im = negate_ymm8r4(xim);
                        return (cv);
                   } 
                   
                   
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cconj_ymm4c8_v2(const ymm4c8_t x) {                                              
                         
                        //register __m256d c;              
                        //c = negate_ymm4c8(*im);
                        //*im = c;
                        ymm4c8_t cv;
                        cv.re = x.re; 
                        cv.im = negate_ymm8r4(x.im);
                        return (cv);
                   } 
                   
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict  csre,
                                       double * __restrict  csim) {

                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_loadu_pd(&xre[0]);
                      ymm1  = _mm256_loadu_pd(&xim[0]);
                      ymm2  = _mm256_mul_pd(_mm256_cos_pd(ymm0),_mm256_cosh_pd(ymm1));
                      _mm256_storeu_pd(&csre[0],ymm2);
                      ymm3  = _mm256_mul_pd(_mm256_sin_pd(ymm0),_mm256_sinh_pd(ymm1));
                      _mm256_storeu_pd(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       double * __restrict  __ATTR_ALIGN__(32) csre,
                                       double * __restrict  __ATTR_ALIGN__(32) csim) {

                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_load_pd(&xre[0]);
                      ymm1  = _mm256_load_pd(&xim[0]);
                      ymm2  = _mm256_mul_pd(_mm256_cos_pd(ymm0),_mm256_cosh_pd(ymm1));
                      _mm256_store_pd(&csre[0],ymm2);
                      ymm3  = _mm256_mul_pd(_mm256_sin_pd(ymm0),_mm256_sinh_pd(ymm1));
                      _mm256_store_pd(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     __m256d * __restrict csre,
                                     __m256d * __restrict csim) {

                      register __m256d ymm0,ymm1;
                      ymm0  = _mm256_mul_pd(_mm256_cos_pd(xre),_mm256_cosh_pd(xim));
                      *csre = ymm0;
                      ymm1  = _mm256_mul_pd(_mm256_sin_pd(xre),_mm256_sinh_pd(xim));
                      *csim = ymm1; 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t ccos_ymm4c8(const __m256d xre,
                                          const __m256d xim) {
                                    
                      ymm4c8_t cv;
                      register __m256d ymm0,ymm1;
                      ymm0  = _mm256_mul_pd(_mm256_cos_pd(xre),_mm256_cosh_pd(xim));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_pd(_mm256_sin_pd(xre),_mm256_sinh_pd(xim));
                      cv.im = ymm1;
                      return (cv); 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t ccos_ymm4c8(const ymm4c8_t x) {
                                    
                      ymm4c8_t cv;
                      register __m256d ymm0,ymm1;
                      ymm0  = _mm256_mul_pd(_mm256_cos_pd(x.re),_mm256_cosh_pd(x.im));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_pd(_mm256_sin_pd(x.re),_mm256_sinh_pd(x.im));
                      cv.im = ymm1;
                      return (cv); 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict  csre,
                                       double * __restrict  csim) {

                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_loadu_pd(&xre[0]);
                      ymm1  = _mm256_loadu_pd(&xim[0]);
                      ymm2  = _mm256_mul_pd(_mm256_cosh_pd(ymm0),_mm256_cos_pd(ymm1));
                      _mm256_storeu_pd(&csre[0],ymm2);
                      ymm3  = _mm256_mul_pd(_mm256_sinh_pd(ymm0),_mm256_sin_pd(ymm1));
                      _mm256_storeu_pd(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       double * __restrict  __ATTR_ALIGN__(32) csre,
                                       double * __restrict  __ATTR_ALIGN__(32) csim) {

                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      ymm0  = _mm256_load_pd(&xre[0]);
                      ymm1  = _mm256_load_pd(&xim[0]);
                      ymm2  = _mm256_mul_pd(_mm256_cosh_pd(ymm0),_mm256_cos_pd(ymm1));
                      _mm256_store_pd(&csre[0],ymm2);
                      ymm3  = _mm256_mul_pd(_mm256_sinh_pd(ymm0),_mm256_sin_pd(ymm1));
                      _mm256_store_pd(&csim[0],ymm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     __m256d * __restrict csre,
                                     __m256d * __restrict csim) {

                      register __m256d ymm0,ymm1;
                      ymm0  = _mm256_mul_pd(_mm256_cosh_pd(xre),_mm256_cos_pd(xim));
                      *csre = ymm0;
                      ymm1  = _mm256_mul_pd(_mm256_sinh_pd(xre),_mm256_sin_pd(xim));
                      *csim = ymm1; 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t ccosh_ymm4c8(const __m256d xre,
                                           const __m256d xim) {
                                          
                      ymm4c8_t cv;
                      register __m256d ymm0,ymm1;
                      ymm0  = _mm256_mul_pd(_mm256_cosh_pd(xre),_mm256_cos_pd(xim));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_pd(_mm256_sinh_pd(xre),_mm256_sin_pd(xim));
                      cv.im = ymm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t ccosh_ymm4c8(const ymm4c8_t x) {
                                          
                      ymm4c8_t cv;
                      register __m256d ymm0,ymm1;
                      ymm0  = _mm256_mul_pd(_mm256_cosh_pd(x.re),_mm256_cos_pd(x.im));
                      cv.re = ymm0;
                      ymm1  = _mm256_mul_pd(_mm256_sinh_pd(x.re),_mm256_sin_pd(x.im));
                      cv.im = ymm1; 
                      return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_ymm4c8(const __m256d xre,
	                             const __m256d xim,
	                             const double n,
	                             __m256d * __restrict powr,
	                             __m256d * __restrict powi) {
	                             
	                register __m256d ymm0,ymm1;
	                register __m256d r,tht;
	                register __m256d vn,pt;
	                register __m256d ta;
	                ymm0  = _mm256_mul_pd(xre,xre);
	                vn    = _mm256_set1_pd(n);
	                ymm1  = _mm256_mul_pd(xim,xim);
	                r     = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
	                tht   = _mm256_atan_pd(_mm256_div_pd(xim,xre));
	                pt    = _mm256_pow_pd(r,vn);
	                ta    = _mm256_mul_pd(vn,tht);
	                *powr = _mm256_mul_pd(pt,_mm256_cos_pd(ta));
	                *powi = _mm256_mul_pd(pt,_mm256_sin_pd(ta));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) pxre,
	                               const double * __restrict __ATTR_ALIGN__(32) pxim,
	                               const double n,
	                               double * __restrict __ATTR_ALIGN__(32) ppowr
	                               double * __restrict __ATTR_ALIGN__(32) ppowi) {
	                  
	                register __m256d xre = _mm256_load_pd(&pxre[0]);
	                register __m256d xim = _mm256_load_pd(&pxim[0]);           
	                register __m256d ymm0,ymm1;
	                register __m256d r,tht;
	                register __m256d vn,pt;
	                register __m256d ta;
	                ymm0  = _mm256_mul_pd(xre,xre);
	                vn    = _mm256_set1_pd(n);
	                ymm1  = _mm256_mul_pd(xim,xim);
	                r     = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
	                tht   = _mm256_atan_pd(_mm256_div_pd(xim,xre));
	                pt    = _mm256_pow_pd(r,vn);
	                ta    = _mm256_mul_pd(vn,tht);
	                _mm256_store_pd(&ppowr[0] ,_mm256_mul_pd(pt,_mm256_cos_pd(ta)));
	                _mm256_store_pd(&ppowi[0] ,_mm256_mul_pd(pt,_mm256_sin_pd(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cpow_ymm4c8_u(const double * __restrict  pxre,
	                               const double * __restrict  pxim,
	                               const double n,
	                               double * __restrict  ppowr
	                               double * __restrict  ppowi) {
	                  
	                register __m256d xre = _mm256_loadu_pd(&pxre[0]);
	                register __m256d xim = _mm256_loadu_pd(&pxim[0]);           
	                register __m256d ymm0,ymm1;
	                register __m256d r,tht;
	                register __m256d vn,pt;
	                register __m256d ta;
	                ymm0  = _mm256_mul_pd(xre,xre);
	                vn    = _mm256_set1_pd(n);
	                ymm1  = _mm256_mul_pd(xim,xim);
	                r     = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
	                tht   = _mm256_atan_pd(_mm256_div_pd(xim,xre));
	                pt    = _mm256_pow_pd(r,vn);
	                ta    = _mm256_mul_pd(vn,tht);
	                _mm256_storeu_pd(&ppowr[0] ,_mm256_mul_pd(pt,_mm256_cos_pd(ta)));
	                _mm256_storeu_pd(&ppowi[0] ,_mm256_mul_pd(pt,_mm256_sin_pd(ta)));                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           ZMM16c4_t cpow_ymm4c8(const ZMM16c4_t x,
	                                  const double n) {
	                   
	                ZMM16c4_t cp;        
	                register __m256d ymm0,ymm1;
	                register __m256d r,tht;
	                register __m256d vn,pt;
	                register __m256d ta;
	                ymm0  = _mm256_mul_pd(x.re,x.re);
	                vn    = _mm256_set1_pd(n);
	                ymm1  = _mm256_mul_pd(x.im,x.im);
	                r     = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
	                tht   = _mm256_atan_pd(_mm256_div_pd(x.im,x.re));
	                pt    = _mm256_pow_pd(r,vn);
	                ta    = _mm256_mul_pd(vn,tht);
	                cp.re = _mm256_mul_pd(pt,_mm256_cos_pd(ta));
	                cp.im = _mm256_mul_pd(pt,_mm256_sin_pd(ta));      
	                return (cp);              
	       }
	       
#include <utility>

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_ymm4c8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim) {
                                      
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      ymm0 = _mm256_loadu_pd(&xre[0]);
                      ymm1 = _mm256_loadu_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_EQ_OQ);
                      ymm2 = _mm256_loadu_pd(&xim[0]);
                      ymm3 = _mm256_loadu_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_EQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                      const double * __restrict __ATTR_ALIGN__(32) xim,
                                      const double * __restrict __ATTR_ALIGN__(32) yre,
                                      const double * __restrict __ATTR_ALIGN__(32) yim) {
                                      
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      ymm0 = _mm256_load_pd(&xre[0]);
                      ymm1 = _mm256_load_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_EQ_OQ);
                      ymm2 = _mm256_load_pd(&xim[0]);
                      ymm3 = _mm256_load_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_EQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_ymm4c8(      const __m256d xre,
                                    const __m256d xim,
                                    const __m256d yre,
                                    const __m256d yim) {
                                    
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm256_cmp_pd_mask(xre,yre,_CMP_EQ_OQ);
                         eqi = _mm256_cmp_pd_mask(xim,yim,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_ymm4c8(      const ymm4c8_t x,
                                    const ymm4c8_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_ymm4c8_u(      const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim) {
                                     
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      ymm0 = _mm256_loadu_pd(&xre[0]);
                      ymm1 = _mm256_loadu_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_GT_OQ);
                      ymm2 = _mm256_loadu_pd(&xim[0]);
                      ymm3 = _mm256_loadu_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_GT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_ymm4c8_a(      const double * __restrict __ATTR_ALIGN__(32) xre,
                                      const double * __restrict __ATTR_ALIGN__(32) xim,
                                      const double * __restrict __ATTR_ALIGN__(32) yre,
                                      const double * __restrict __ATTR_ALIGN__(32) yim) {
                                      
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      ymm0 = _mm256_load_pd(&xre[0]);
                      ymm1 = _mm256_load_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_GT_OQ);
                      ymm2 = _mm256_load_pd(&xim[0]);
                      ymm3 = _mm256_load_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_GT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_ymm4c8(      const __m256d xre,
                                    const __m256d xim,
                                    const __m256d yre,
                                    const __m256d yim) {
                            
                         __mmask8 eqr;
                         __mmask8 eqi;      
                         eqr = _mm256_cmp_pd_mask(xre,yre,_CMP_GT_OQ);
                         eqi = _mm256_cmp_pd_mask(xim,yim,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_ymm4c8(      const ymm4c8_t x,
                                    const ymm4c8_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_ymm4c8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim){
                                      
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      ymm0 = _mm256_loadu_pd(&xre[0]);
                      ymm1 = _mm256_loadu_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_LT_OQ);
                      ymm2 = _mm256_loadu_pd(&xim[0]);
                      ymm3 = _mm256_loadu_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_LT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_ymm4c8_a(      const double * __restrict __ATTR_ALIGN__(32) xre,
                                      const double * __restrict __ATTR_ALIGN__(32) xim,
                                      const double * __restrict __ATTR_ALIGN__(32) yre,
                                      const double * __restrict __ATTR_ALIGN__(32) yim) {
                                     
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      ymm0 = _mm256_load_pd(&xre[0]);
                      ymm1 = _mm256_load_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_LT_OQ);
                      ymm2 = _mm256_load_pd(&xim[0]);
                      ymm3 = _mm256_load_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_LT_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_ymm4c8(      const __m256d xre,
                                    const __m256d xim,
                                    const __m256d yre,
                                    const __m256d yim) {
                               
                         __mmask8 eqr;
                         __mmask8 eqi;      
                         eqr = _mm256_cmp_pd_mask(xre,yre,_CMP_LT_OQ);
                         eqi = _mm256_cmp_pd_mask(xim,yim,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_ymm4c8(     const ymm4c8_t x,
                                    const ymm4c8_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_ymm4c8_u(     const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim) {
                                      
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi; 
                      ymm0 = _mm256_loadu_pd(&xre[0]);
                      ymm1 = _mm256_loadu_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_NEQ_OQ);
                      ymm2 = _mm256_loadu_pd(&xim[0]);
                      ymm3 = _mm256_loadu_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_NEQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                      const double * __restrict __ATTR_ALIGN__(32) xim,
                                      const double * __restrict __ATTR_ALIGN__(32) yre,
                                      const double * __restrict __ATTR_ALIGN__(32) yim) {
                                      
                      register __m256d ymm0,ymm1,ymm2,ymm3;
                      __mmask8 eqr;
                      __mmask8 eqi;
                      ymm0 = _mm256_load_pd(&xre[0]);
                      ymm1 = _mm256_load_pd(&yre[0]);
                      eqr  = _mm256_cmp_pd_mask(ymm0,ymm1,_CMP_NEQ_OQ);
                      ymm2 = _mm256_load_pd(&xim[0]);
                      ymm3 = _mm256_load_pd(&yim[0]);
                      eqi  = _mm256_cmp_pd_mask(ymm2,ymm3,_CMP_NEQ_OQ);
                      return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_ymm4c8(     const __m256d xre,
                                    const __m256d xim,
                                    const __m256d yre,
                                    const __m256d yim) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm256_cmp_pd_mask(xre,yre,_CMP_NEQ_OQ);
                         eqi = _mm256_cmp_pd_mask(xim,yim,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_ymm4c8(      const ymm4c8_t x,
                                     const ymm4c8_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict cexpr,
                                       double * __restrict cexpi ) {

                        register const __m256d I = _mm256_set1_pd(1.0f);
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_loadu_pd(&xre[0]);
                        ymm1  = _mm256_loadu_pd(&xim[0]);
                        ymm2  = _mm256_exp_pd(ymm0);
                        ymm3  = _mm256_mul_pd(ymm2,_mm256_cos_pd(ymm1));
                        _mm256_storeu_pd(&cexpr[0],ymm3);
                        ymm4  = _mm256_mul_pd(ymm2,_mm256_mul_pd(_mm256_sin_pd(ymm1),I));
                        _mm256_storeu_pd(&cexpi[0],ymm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       double * __restrict __ATTR_ALIGN__(32) cexpr,
                                       double * __restrict __ATTR_ALIGN__(32) cexpi ) {

                        register const __m256d I = _mm256_set1_pd(1.0f);
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        ymm0  = _mm256_load_pd(&xre[0]);
                        ymm1  = _mm256_load_pd(&xim[0]);
                        ymm2  = _mm256_exp_pd(ymm0);
                        ymm3  = _mm256_mul_pd(ymm2,_mm256_cos_pd(ymm1));
                        _mm256_store_pd(&cexpr[0],ymm3);
                        ymm4  = _mm256_mul_pd(ymm2,_mm256_mul_pd(_mm256_sin_pd(ymm1),I));
                        _mm256_store_pd(&cexpi[0],ymm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_ymm4c8(const __m256d xre,
                                     const __m256d xim,
                                     __m256d * __restrict cexpr,
                                     __m256d * __restrict cexpi) {

                        register const __m256d I = _mm256_set1_pd(1.0f);
                        register __m256d ymm0;
                        ymm0   = _mm256_exp_pd(xre);
                        *cexpr = _mm256_mul_pd(ymm0,_mm256_cos_pd(xim));
                        *cexpi = _mm256_mul_pd(ymm0,_mm256_mul_pd(_mm256_sin_pd(xim),I));
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cexp_ymm4c8(const __m256d xre,
                                          const __m256d xim) {
                                     
                        ymm4c8_t cv;
                        register const __m256d I = _mm256_set1_pd(1.0f);
                        register __m256d ymm0;
                        ymm0   = _mm256_exp_pd(xre);
                        cv.re = _mm256_mul_pd(ymm0,_mm256_cos_pd(xim));
                        cv.im = _mm256_mul_pd(ymm0,_mm256_mul_pd(_mm256_sin_pd(xim),I));
                        return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cexp_ymm4c8(const ymm4c8_t x) {
                                     
                        ymm4c8_t cv;
                        register const __m256d I = _mm256_set1_pd(1.0f);
                        register __m256d ymm0;
                        ymm0   = _mm256_exp_pd(x.re);
                        cv.re = _mm256_mul_pd(ymm0,_mm256_cos_pd(x.im));
                        cv.im = _mm256_mul_pd(ymm0,_mm256_mul_pd(_mm256_sin_pd(x.im),I));
                        return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_ymm4c8_u(const double * __restrict rho,
                                         const double * __restrict tht,
                                         double * __restrict  re,
                                         double * __restrict  im) {

                         register __m256d ymm0,ymm1,ymm2,ymm3;
                         ymm0 = _mm256_loadu_pd(&rho[0]);
                         ymm1 = _mm256_loadu_pd(&tht[0]);
                         ymm2 = _mm256_mul_pd(ymm0,_mm256_cos_pd(ymm1)); //tht
                         _mm256_storeu_pd(&re[0],ymm2);
                         ymm3 = _mm256_mul_pd(ymm0,_mm256_sin_pd(ymm1)); //tht
                         _mm256_storeu_pd(&im[0],ymm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) rho,
                                         const double * __restrict __ATTR_ALIGN__(32) tht,
                                         double * __restrict  __ATTR_ALIGN__(32) re,
                                         double * __restrict  __ATTR_ALIGN__(32) im) {

                         register __m256d ymm0,ymm1,ymm2,ymm3;
                         ymm0 = _mm256_load_pd(&rho[0]);
                         ymm1 = _mm256_load_pd(&tht[0]);
                         ymm2 = _mm256_mul_pd(ymm0,_mm256_cos_pd(ymm1)); //tht
                         _mm256_store_pd(&re[0],ymm2);
                         ymm3 = _mm256_mul_pd(ymm0,_mm256_sin_pd(ymm1)); //tht
                         _mm256_store_pd(&im[0],ymm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_ymm4c8(const __m256d rho,
                                       const __m256d tht,
                                       __m256d * __restrict re,
                                       __m256d * __restrict im) {

                        register __m256d ymm0,ymm1;
                        ymm0 = _mm256_mul_pd(rho,_mm256_cos_pd(tht));
                        *re  = ymm0;
                        ymm1 = _mm256_mul_pd(rho,_mm256_sin_pd(tht));
                        *im  = ymm1;
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cpolar_ymm4c8(const __m256d rho,
                                            const __m256d tht) {
                                      
                        ymm4c8_t cv
                        register __m256d ymm0,ymm1;
                        ymm0 = _mm256_mul_pd(rho,_mm256_cos_pd(tht));
                        cv.re  = ymm0;
                        ymm1 = _mm256_mul_pd(rho,_mm256_sin_pd(tht));
                        cv.im  = ymm1;
                        return (cv);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_ymm4c8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict wrkc,
                                       double * __restrict csqr,
                                       double * __restrict csqi) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d half = _mm256_set1_pd(0.5f);
                        cabs_ymm4c8_u(xre,xim,wrkc);
                        ymm0  = _mm256_loadu_pd(&xre[0]);
                        ymm1  = _mm256_loadu_pd(&wrkc[0]);
                        ymm2  = _mm256_mul_pd(half,_mm256_add_pd(ymm1,ymm0));
                        _mm256_storeu_pd(&csqr[0],_mm256_sqrt_pd(ymm2));
                        ymm3  = _mm256_mul_pd(half,_mm256_sub_pd(ymm1,ymm0));
                        _mm256_storeu_pd(&csqi[0],_mm256_sqrt_pd(ymm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                       const double * __restrict __ATTR_ALIGN__(32) xim,
                                       double * __restrict __ATTR_ALIGN__(32) wrkc,
                                       double * __restrict __ATTR_ALIGN__(32) csqr,
                                       double * __restrict __ATTR_ALIGN__(32) csqi) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d half = _mm256_set1_pd(0.5f);
                        cabs_ymm4c8_a(xre,xim,wrkc);
                        ymm0  = _mm256_load_pd(&xre[0]);
                        ymm1  = _mm256_load_pd(&wrkc[0]);
                        ymm2  = _mm256_mul_pd(half,_mm256_add_pd(ymm1,ymm0));
                        _mm256_store_pd(&csqr[0],_mm256_sqrt_pd(ymm2));
                        ymm3  = _mm256_mul_pd(half,_mm256_sub_pd(ymm1,ymm0));
                        _mm256_store_pd(&csqi[0],_mm256_sqrt_pd(ymm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_ymm4c8(const __m256d xre,
                                      const __m256d xim,
                                      __m256d * __restrict wrkc,
                                      __m256d * __restrict csqr,
                                      __m256d * __restrict csqi) {

                       register __m256d ymm0,ymm1;
                       register __m256d half = _mm256_set1_pd(0.5f); 
                       cabs_ymm4c8(xre,xim,wrkc);
                       ymm0  = _mm256_mul_pd(half,_mm256_add_pd(*wrkc,xre));
                       *csqr = ymm0;
                       ymm1  = _mm256_mul_pd(half,_mm256_sub_pd(*wrkc,xre));
                       *csqi = ymm1; 
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t csqrt_ymm4c8(const __m256d xre,
                                           const __m256d xim,
                                          __m256d * __restrict wrkc) {
                                          
                       ymm4c8_t cv;
                       register __m256d ymm0,ymm1;
                       register __m256d half = _mm256_set1_pd(0.5f); 
                       cabs_ymm4c8(xre,xim,wrkc);
                       ymm0  = _mm256_mul_pd(half,_mm256_add_pd(*wrkc,xre));
                       cv.re = ymm0;
                       ymm1  = _mm256_mul_pd(half,_mm256_sub_pd(*wrkc,xre));
                       cv.im = ymm1; 
                       return (cv);
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t csqrt_ymm4c8(const ymm4c8_t x,
                                          __m256d * __restrict wrkc) {
                                          
                       ymm4c8_t cv;
                       register __m256d ymm0,ymm1;
                       register __m256d half = _mm256_set1_pd(0.5f); 
                       cabs_ymm4c8(x.re,x.im,wrkc);
                       ymm0  = _mm256_mul_pd(half,_mm256_add_pd(*wrkc,x.re));
                       cv.re = ymm0;
                       ymm1  = _mm256_mul_pd(half,_mm256_sub_pd(*wrkc,x.re));
                       cv.im = ymm1; 
                       return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_ymm4c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict zre,
                                             double * __restrict zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        rep  = _mm256_fmsub_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        imp  = _mm256_fmadd_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        ymm0 = _mm256_mul_pd(rep,rep);
                        ymm1 = _mm256_mul_pd(imp,imp);
                        ymm2 = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        _mm256_storeu_pd(&zre[0], _mm256_div_pd(rep,ymm2));
                        _mm256_storeu_pd(&zim[0], _mm256_div_pd(imp,ymm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                             const double * __restrict __ATTR_ALIGN__(32) xim,
                                             const double * __restrict __ATTR_ALIGN__(32) yre,
                                             const double * __restrict __ATTR_ALIGN__(32) yim,
                                             double * __restrict __ATTR_ALIGN__(32) zre,
                                             double * __restrict __ATTR_ALIGN__(32) zim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        rep  = _mm256_fmsub_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        imp  = _mm256_fmadd_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        ymm0 = _mm256_mul_pd(rep,rep);
                        ymm1 = _mm256_mul_pd(imp,imp);
                        ymm2 = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        _mm256_store_pd(&zre[0], _mm256_div_pd(rep,ymm2));
                        _mm256_store_pd(&zim[0], _mm256_div_pd(imp,ymm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_ymm4c8(  const __m256d  xre,
                                             const __m256d  xim,
                                             const __m256d  yre,
                                             const __m256d  yim,
                                             __m256d * __restrict zre,
                                             __m256d * __restrict zim) {

                        register __m256d rep,imp,ymm0,ymm1,ymm2;
                        rep  = _mm256_fmsub_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        imp  = _mm256_fmadd_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        ymm0 = _mm256_mul_pd(rep,rep);
                        ymm1 = _mm256_mul_pd(imp,imp);
                        ymm2 = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        *zre = _mm256_div_pd(rep,ymm2);
                        *zim = _mm256_div_pd(imp,ymm2);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cnorm_prod_ymm4c8(  const __m256d  xre,
                                                  const __m256d  xim,
                                                  const __m256d  yre,
                                                  const __m256d  yim) {
                                             
                        ymm4c8_t cv;
                        register __m256d rep,imp,ymm0,ymm1,ymm2;
                        rep  = _mm256_fmsub_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        imp  = _mm256_fmadd_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        ymm0 = _mm256_mul_pd(rep,rep);
                        ymm1 = _mm256_mul_pd(imp,imp);
                        ymm2 = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        cv.re = _mm256_div_pd(rep,ymm2);
                        cv.im = _mm256_div_pd(imp,ymm2);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cnorm_prod_ymm4c8(  const ymm4c8_t x,
                                                  const ymm4c8_t y) {
                                             
                        ymm4c8_t cv;
                        register __m256d rep,imp,ymm0,ymm1,ymm2;
                        rep  = _mm256_fmsub_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        imp  = _mm256_fmadd_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        ymm0 = _mm256_mul_pd(rep,rep);
                        ymm1 = _mm256_mul_pd(imp,imp);
                        ymm2 = _mm256_sqrt_pd(_mm256_add_pd(ymm0,ymm1));
                        cv.re = _mm256_div_pd(rep,ymm2);
                        cv.im = _mm256_div_pd(imp,ymm2);
                        return (cv);
             }
             
#include "GMS_simd_utils.hpp"


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm4c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm256_fmsub_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                             const double * __restrict __ATTR_ALIGN__(32) xim,
                                             const double * __restrict __ATTR_ALIGN__(32) yre,
                                             const double * __restrict __ATTR_ALIGN__(32) yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm256_fmsub_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
              } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm4c8(const __m256d xre,
                                           const __m256d xim,
                                           const __m256d yre,
                                           const __m256d yim,
                                           double * __restrict mre,
                                           double * __restrict min) {

                        register __m256d rep,imp;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        sre = 0.0f;
                        rep  = _mm256_fmsub_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_ymm4c8(const ymm4c8_t x,
                                           const ymm4c8_t y,
                                           double * __restrict mre,
                                           double * __restrict min) {

                        register __m256d rep,imp;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        sre = 0.0f;
                        rep  = _mm256_fmsub_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        sre  = ymm8r4_horizontal_sum(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm256_fmadd_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        sim  = ymm8r4_horizontal_sum(imp);
                        *mim = sim*inv16;
             }


               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm4c8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        imp  = _mm256_fmadd_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_pd(ymm1,ymm1,
                                               _mm256_mul_pd(ymm3,ymm3));
                        rquot = _mm256_div_pd(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_pd(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                             const double * __restrict __ATTR_ALIGN__(32) xim,
                                             const double * __restrict __ATTR_ALIGN__(32) yre,
                                             const double * __restrict __ATTR_ALIGN__(32) yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        imp  = _mm256_fmadd_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_pd(ymm1,ymm1,
                                               _mm256_mul_pd(ymm3,ymm3));
                        rquot = _mm256_div_pd(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_pd(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm4c8(  const __m256d xre,
                                             const __m256d xim,
                                             const __m256d yre,
                                             const __m256d yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        imp  = _mm256_fmadd_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_pd(yre,yre,
                                               _mm256_mul_pd(yim,yim));
                        rquot = _mm256_div_pd(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_pd(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_ymm4c8(  const ymm4c8_t x,
                                             const ymm4c8_t y,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        sre  = 0.0f;
                        rep  = _mm256_fmsub_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        imp  = _mm256_fmadd_pd(xim,yre,
                                               _mm256_mul_pd(x.re,y.im));
                        sim  = 0.0f;
                        den  = _mm256_fmadd_pd(y.re,y.re,
                                               _mm256_mul_pd(y.im,y.im));
                        rquot = _mm256_div_pd(rep,den);
                        sre   = ymm8r4_horizontal_sum(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm256_div_pd(imp,den);
                        sim   = ymm8r4_horizontal_sum(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm4c8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim ) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp,magc1,magc2,vcmag;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        rep  = _mm256_fmadd_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        _mm256_storeu_pd(&mre[0], _mm256_div_pd(rep,vcmag));
                        _mm256_storeu_pd(&mim[0], _mm256_div_pd(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                              const double * __restrict __ATTR_ALIGN__(32) xim,
                                              const double * __restrict __ATTR_ALIGN__(32) yre,
                                              const double * __restrict __ATTR_ALIGN__(32) yim,
                                              double * __restrict __ATTR_ALIGN__(32) mre,
                                              double * __restrict __ATTR_ALIGN__(32) mim ) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d rep,imp,magc1,magc2,vcmag;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        rep  = _mm256_fmad_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        _mm256_store_pd(&mre[0], _mm256_div_pd(rep,vcmag));
                        _mm256_store_pd(&mim[0], _mm256_div_pd(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm4c8(const __m256d xre,
                                            const __m256d xim,
                                            const __m256d yre,
                                            const __m256d yim,
                                            __m256d * __restrict mre,
                                            __m256d * __restrict mim) {

                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        *mre = _mm256_div_pd(rep,vcmag);
                        *mim = _mm256_div_pd(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_ymm4c8(const ymm4c8_t x,
                                            const ymm4c8_t y,
                                            __m256d * __restrict mre,
                                            __m256d * __restrict mim) {

                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        *mre = _mm256_div_pd(rep,vcmag);
                        *mim = _mm256_div_pd(imp,vcmag);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cnorm_cprod_ymm4c8(const __m256d xre,
                                                 const __m256d xim,
                                                 const __m256d yre,
                                                 const __m256d yim) {
                                               
                        ymm4c8_t cv;
                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        cv.re = _mm256_div_pd(rep,vcmag);
                        cv.im = _mm256_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cnorm_cprod_ymm4c8(const ymm4c8_t x,
                                                 const ymm4c8_t y) {
                                               
                        ymm4c8_t cv;
                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmad_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        cv.re = _mm256_div_pd(rep,vcmag);
                        cv.im = _mm256_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
             


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm4c8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim) {
                      
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d re,im;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        re   = _mm256_fmadd_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                              const double * __restrict __ATTR_ALIGN__(32) xim,
                                              const double * __restrict __ATTR_ALIGN__(32) yre,
                                              const double * __restrict __ATTR_ALIGN__(32) yim,
                                              double * __restrict mre,
                                              double * __restrict mim) {
                      
                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d re,im;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        re   = _mm256_fmadd_pd(ymm0,ymm1,
                                               _mm256_mul_pd(ymm2,ymm3));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_pd(ymm2,ymm1,
                                               _mm256_mul_pd(ymm0,ymm3));
                        sim  = ymm8r4_horizontal_sum(rep);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm4c8(const __m256d xre,
                                            const __m256d xim,
                                            const __m256d yre,
                                            const __m256d yim,
                                            double * __restrict mre,
                                            double * __restrict min) {

                        register __m256d re,im;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        re   = _mm256_fmadd_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_ymm4c8(const ymm4c8_t x,
                                            const ymm4c8_t y,
                                            double * __restrict mre,
                                            double * __restrict min) {

                        register __m256d re,im;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        re   = _mm256_fmadd_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_fmsub_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm4c8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        register __m256d re,im;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        re   = _mm256_loadu_pd(&xre[0]);
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_loadu_pd(&xim[0]);
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                              const double * __restrict __ATTR_ALIGN__(32) xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        register __m256d re,im;
                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
                        re   = _mm256_load_pd(&xre[0]);
                        sre  = ymm8r4_horizontal_sum(re);
                        *mre = sre*inv16;
                        im   = _mm256_load_pd(&xim[0]);
                        sim  = ymm8r4_horizontal_sum(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_ymm4c8(  const __m256d xre,
                                              const __m256d xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
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
                   void arith_cmean_ymm4c8(  const ymm4c8_t x,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        constexpr double inv16 = 0.0625f;
                        double sre,sim;
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
                   void cnormalize_ymm4c8_u( const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim ) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d re,im,cvmag;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(ymm0,ymm1,
                                                              _mm256_mul_pd(ymm2,ymm3)));
                        _mm256_storeu_pd(&mre[0], _mm256_div_pd(ymm0,cvmag));
                        _mm256_storeu_pd(&mim[0], _mm256_div_pd(ymm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm4c8_a( const double * __restrict __ATTR_ALIGN__(32) xre,
                                              const double * __restrict __ATTR_ALIGN__(32) xim,
                                              const double * __restrict __ATTR_ALIGN__(32) yre,
                                              const double * __restrict __ATTR_ALIGN__(32) yim,
                                              double * __restrict __ATTR_ALIGN__(32) mre,
                                              double * __restrict __ATTR_ALIGN__(32) mim ) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d re,im,cvmag;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(ymm0,ymm1,
                                                              _mm256_mul_pd(ymm2,ymm3)));
                        _mm256_store_pd(&mre[0], _mm256_div_pd(ymm0,cvmag));
                        _mm256_store_pd(&mim[0], _mm256_div_pd(ymm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm4c8( const __m256d xre,
                                            const __m256d xim,
                                            const __m256d yre,
                                            const __m256d yim,
                                            __m256d * __restrict mre,
                                            __m256d * __restrict mim ) {

                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(xre,yre,
                                                    _mm256_mul_pd(xim,yim)));
                        *mre = _mm256_div_pd(xre,cvmag));
                        *mim =  _mm256_div_pd(xim,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_ymm4c8( const ymm4c8_t x,
                                            const ymm4c8_t y,
                                            __m256d * __restrict mre,
                                            __m256d * __restrict mim ) {

                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(x.re,y.re,
                                                    _mm256_mul_pd(x.im,y.im)));
                        *mre = _mm256_div_pd(x.re,cvmag));
                        *mim =  _mm256_div_pd(x.im,cvmag));
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cnormalize_ymm4c8( const __m256d xre,
                                                 const __m256d xim,
                                                 const __m256d yre,
                                                 const __m256d yim) {
                                            
                        ymm4c8_t cv;
                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(xre,yre,
                                                    _mm256_mul_pd(xim,yim)));
                        cv.re = _mm256_div_pd(xre,cvmag));
                        cv.im =  _mm256_div_pd(xim,cvmag));
                        return (cv);
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   ymm4c8_t cnormalize_ymm4c8( const ymm4c8_t x,
                                                 const ymm4c8_t y,) {
                                            
                        ymm4c8_t cv;
                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(x.re,y.re,
                                                    _mm256_mul_pd(x.im,y.im)));
                        cv.re = _mm256_div_pd(x.re,cvmag));
                        cv.im =  _mm256_div_pd(x.im,cvmag));
                        return (cv);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm4c8_u( const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d cvmag;
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&yre[0]);
                        ymm2 = _mm256_loadu_pd(&xim[0]);
                        ymm3 = _mm256_loadu_pd(&yim[0]);
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(ymm0,ymm1,
                                                          _mm256_mul_pd(ymm2,ymm3)));
                        _mm256_storeu_pd(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm4c8_a( const double * __restrict __ATTR_ALIGN__(32) xre,
                                              const double * __restrict __ATTR_ALIGN__(32) xim,
                                              const double * __restrict __ATTR_ALIGN__(32) yre,
                                              const double * __restrict __ATTR_ALIGN__(32) yim,
                                              double * __restrict __ATTR_ALIGN__(32) mre) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d cvmag;
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&yre[0]);
                        ymm2 = _mm256_load_pd(&xim[0]);
                        ymm3 = _mm256_load_pd(&yim[0]);
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(ymm0,ymm1,
                                                          _mm256_mul_pd(ymm2,ymm3)));
                        _mm256_store_pd(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm4c8(   const __m256d xre,
                                              const __m256d xim,
                                              const __m256d yre,
                                              const __m256d yim,
                                              __m256d * __restrict  mre) {

                        register __m256d cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(xre,yre,
                                                          _mm256_mul_pd(xim,yim)));
                        *mre = cvmag;
             }
             
             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_ymm4c8(   const ymm4c8_t x,
                                              const ymm4c8_t y,
                                              __m256d * __restrict  mre) {

                        register __m256d cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(x.re,y.re,
                                                          _mm256_mul_pd(x.im,y.im)));
                        *mre = cvmag;
             }


#include <complex>


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_2xr4_c4_unroll16x(double * __restrict __ATTR_ALIGN__(32) xre,
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
                   void copy_c4_2xr4_unroll16x( double * __restrict __ATTR_ALIGN__(32) xre,
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
                   void copy_2xr4_c4_ymm4c8_u(const double * __restrict xre,
                                               const double * __restrict xim,
                                               std::complex<double> * __restrict yc) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d ymm4,ymm5;
                        double * __restrict pyc = reinterpret_cast<double*>(&yc[0]);
                        ymm0 = _mm256_loadu_pd(&xre[0]);
                        ymm1 = _mm256_loadu_pd(&xim[0]);
                        ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                        ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                        ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                        _mm256_storeu_pd(&pyc[0],ymm4);
                        ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                        _mm256_storeu_pd(&pyc[8],ymm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_ymm4c8(  const __m256d xre,
                                               const __m256d xim,
                                               std::complex<double> * __restrict yc) {

                        register __m256d ymm2,ymm3;
                        register __m256d ymm4,ymm5;
                        double * __restrict pyc = reinterpret_cast<double*>(&yc[0]);
                        ymm2 = _mm256_unpacklo_pd(xre,xim);
                        ymm3 = _mm256_unpackhi_pd(xre,xim);
                        ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                        _mm256_storeu_pd(&pyc[0],ymm4);
                        ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                        _mm256_storeu_pd(&pyc[8],ymm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_ymm4c8_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                               const double * __restrict __ATTR_ALIGN__(32) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(32) yc) {

                        register __m256d ymm0,ymm1,ymm2,ymm3;
                        register __m256d ymm4,ymm5;
                        double * __restrict pyc = reinterpret_cast<double*>(&yc[0]);
                        ymm0 = _mm256_load_pd(&xre[0]);
                        ymm1 = _mm256_load_pd(&xim[0]);
                        ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                        ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                        ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                        _mm256_store_pd(&pyc[0],ymm4);
                        ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                        _mm256_store_pd(&pyc[8],ymm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void copy_2r4c4_ymm4c8_blocked_a(const double * __restrict __ATTR_ALIGN__(32) xre,
                                                   const double * __restrict __ATTR_ALIGN__(32) xim,
                                                   std::complex<double> * __restrict __ATTR_ALIGN__(32) yc,
                                                   const int32_t n) { // size of array std::complex<double> len elements.

                        double * __restrict pyc = reinterpret_cast<double*>(&yc[0]);
                        if(n <= 1) {
                            return;
                        }
                        else if(n <= 8) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_load_pd(&xre[0]);
                            ymm1 = _mm256_load_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[8],ymm5);
                        }
                        else if(n <= 16) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_load_pd(&xre[0]);
                            ymm1 = _mm256_load_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[8]);
                            ymm1 = _mm256_load_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[24],ymm5);
                       }
                       else if(n <= 32) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_load_pd(&xre[0]);
                            ymm1 = _mm256_load_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[8]);
                            ymm1 = _mm256_load_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[16]);
                            ymm1 = _mm256_load_pd(&xim[16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[24]);
                            ymm1 = _mm256_load_pd(&xim[24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[56],ymm5);
                       }
                       else if(n <= 64) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_load_pd(&xre[0]);
                            ymm1 = _mm256_load_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[8]);
                            ymm1 = _mm256_load_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[16]);
                            ymm1 = _mm256_load_pd(&xim[16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[24]);
                            ymm1 = _mm256_load_pd(&xim[24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[56],ymm5);
                            ymm0 = _mm256_load_pd(&xre[32]);
                            ymm1 = _mm256_load_pd(&xim[32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[72],ymm5);
                            ymm0 = _mm256_load_pd(&xre[40]);
                            ymm1 = _mm256_load_pd(&xim[40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[88],ymm5);
                            ymm0 = _mm256_load_pd(&xre[48]);
                            ymm1 = _mm256_load_pd(&xim[48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[104],ymm5);
                            ymm0 = _mm256_load_pd(&xre[56]);
                            ymm1 = _mm256_load_pd(&xim[56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[120],ymm5);
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
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_load_pd(&xre[0]);
                            ymm1 = _mm256_load_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[8]);
                            ymm1 = _mm256_load_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[16]);
                            ymm1 = _mm256_load_pd(&xim[16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[24]);
                            ymm1 = _mm256_load_pd(&xim[24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[56],ymm5);
                            ymm0 = _mm256_load_pd(&xre[32]);
                            ymm1 = _mm256_load_pd(&xim[32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[72],ymm5);
                            ymm0 = _mm256_load_pd(&xre[40]);
                            ymm1 = _mm256_load_pd(&xim[40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[88],ymm5);
                            ymm0 = _mm256_load_pd(&xre[48]);
                            ymm1 = _mm256_load_pd(&xim[48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[104],ymm5);
                            ymm0 = _mm256_load_pd(&xre[56]);
                            ymm1 = _mm256_load_pd(&xim[56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[120],ymm5);
                            ymm0 = _mm256_load_pd(&xre[64]);
                            ymm1 = _mm256_load_pd(&xim[64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[136],ymm5);
                            ymm0 = _mm256_load_pd(&xre[72]);
                            ymm1 = _mm256_load_pd(&xim[72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[152],ymm5);
                            ymm0 = _mm256_load_pd(&xre[80]);
                            ymm1 = _mm256_load_pd(&xim[80]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[168],ymm5);
                            ymm0 = _mm256_load_pd(&xre[88]);
                            ymm1 = _mm256_load_pd(&xim[88]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[184],ymm5);
                            ymm0 = _mm256_load_pd(&xre[96]);
                            ymm1 = _mm256_load_pd(&xim[96]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[200],ymm5);
                            ymm0 = _mm256_load_pd(&xre[104]);
                            ymm1 = _mm256_load_pd(&xim[104]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[216],ymm5);
                            ymm0 = _mm256_load_pd(&xre[112]);
                            ymm1 = _mm256_load_pd(&xim[112]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[232],ymm5);
                            ymm0 = _mm256_load_pd(&xre[120]);
                            ymm1 = _mm256_load_pd(&xim[120]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[248],ymm5);
                       }
                      else if(n > 128) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
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
                            ymm0 = _mm256_load_pd(&xre[i+0]);
                            ymm1 = _mm256_load_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+8]);
                            ymm1 = _mm256_load_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+16]);
                            ymm1 = _mm256_load_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+24]);
                            ymm1 = _mm256_load_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+56],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+32]);
                            ymm1 = _mm256_load_pd(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+72],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+40]);
                            ymm1 = _mm256_load_pd(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+88],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+48]);
                            ymm1 = _mm256_load_pd(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+104],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+56]);
                            ymm1 = _mm256_load_pd(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+120],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+64]);
                            ymm1 = _mm256_load_pd(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+136],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+72]);
                            ymm1 = _mm256_load_pd(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+152],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+80]);
                            ymm1 = _mm256_load_pd(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+168],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+88]);
                            ymm1 = _mm256_load_pd(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+184],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+96]);
                            ymm1 = _mm256_load_pd(&xim[i+96]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+200],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+104]);
                            ymm1 = _mm256_load_pd(&xim[i+104]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+216],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+112]);
                            ymm1 = _mm256_load_pd(&xim[i+112]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+232],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+120]);
                            ymm1 = _mm256_load_pd(&xim[i+120]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+248],ymm5);
                     }

                       for(; (i+95) < n; i += 96) {
                            ymm0 = _mm256_load_pd(&xre[i+0]);
                            ymm1 = _mm256_load_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+8]);
                            ymm1 = _mm256_load_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+16]);
                            ymm1 = _mm256_load_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+24]);
                            ymm1 = _mm256_load_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+56],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+32]);
                            ymm1 = _mm256_load_pd(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+72],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+40]);
                            ymm1 = _mm256_load_pd(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+88],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+48]);
                            ymm1 = _mm256_load_pd(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+104],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+56]);
                            ymm1 = _mm256_load_pd(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+120],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+64]);
                            ymm1 = _mm256_load_pd(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+136],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+72]);
                            ymm1 = _mm256_load_pd(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+152],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+80]);
                            ymm1 = _mm256_load_pd(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+168],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+88]);
                            ymm1 = _mm256_load_pd(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+184],ymm5);
                     }    

                       for(; (i+79) < n; i += 80) {
                            ymm0 = _mm256_load_pd(&xre[i+0]);
                            ymm1 = _mm256_load_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+8]);
                            ymm1 = _mm256_load_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+16]);
                            ymm1 = _mm256_load_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+24]);
                            ymm1 = _mm256_load_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+56],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+32]);
                            ymm1 = _mm256_load_pd(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+72],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+40]);
                            ymm1 = _mm256_load_pd(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+88],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+48]);
                            ymm1 = _mm256_load_pd(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+104],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+56]);
                            ymm1 = _mm256_load_pd(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+120],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+64]);
                            ymm1 = _mm256_load_pd(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+136],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+72]);
                            ymm1 = _mm256_load_pd(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+152],ymm5);
                     }

                       for(; (i+31) < n; i += 32) {
                            ymm0 = _mm256_load_pd(&xre[i+0]);
                            ymm1 = _mm256_load_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+8]);
                            ymm1 = _mm256_load_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+16]);
                            ymm1 = _mm256_load_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_load_pd(&xre[i+24]);
                            ymm1 = _mm256_load_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+56],ymm5);
                     }

                     for(; (i+15) < n; i += 16) {
                            ymm0 = _mm256_load_pd(&xre[i+0]);
                            ymm1 = _mm256_load_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+8],ymm5);
                            
                     }

                       for(; (i+0) < n; i += 1) {
                             const double re = xre[i];
                             const double im = xim[i];
                             pyc[i]         = re;
                             pyc[i+1]       = im;
                      
                       }

                }


 
            }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void copy_2r4c4_ymm4c8_blocked_u(const double * __restrict  xre,
                                                   const double * __restrict  xim,
                                                   std::complex<double> *  yc,
                                                   const int32_t n) { // size of array std::complex<double> len elements.

                        double * __restrict pyc = reinterpret_cast<double*>(&yc[0]);
                        if(n <= 1) {
                            return;
                        }
                        else if(n <= 8) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_loadu_pd(&xre[0]);
                            ymm1 = _mm256_loadu_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[8],ymm5);
                        }
                        else if(n <= 16) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_loadu_pd(&xre[0]);
                            ymm1 = _mm256_loadu_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[8]);
                            ymm1 = _mm256_loadu_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[24],ymm5);
                       }
                       else if(n <= 32) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_loadu_pd(&xre[0]);
                            ymm1 = _mm256_loadu_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[8]);
                            ymm1 = _mm256_loadu_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[16]);
                            ymm1 = _mm256_loadu_pd(&xim[16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[24]);
                            ymm1 = _mm256_loadu_pd(&xim[24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[56],ymm5);
                       }
                       else if(n <= 64) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_loadu_pd(&xre[0]);
                            ymm1 = _mm256_loadu_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[8]);
                            ymm1 = _mm256_loadu_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[16]);
                            ymm1 = _mm256_loadu_pd(&xim[16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[24]);
                            ymm1 = _mm256_loadu_pd(&xim[24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[56],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[32]);
                            ymm1 = _mm256_loadu_pd(&xim[32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[72],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[40]);
                            ymm1 = _mm256_loadu_pd(&xim[40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[88],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[48]);
                            ymm1 = _mm256_loadu_pd(&xim[48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[104],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[56]);
                            ymm1 = _mm256_loadu_pd(&xim[56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[120],ymm5);
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
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
                            ymm0 = _mm256_loadu_pd(&xre[0]);
                            ymm1 = _mm256_loadu_pd(&xim[0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[8]);
                            ymm1 = _mm256_loadu_pd(&xim[8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[16]);
                            ymm1 = _mm256_loadu_pd(&xim[16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[24]);
                            ymm1 = _mm256_loadu_pd(&xim[24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[56],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[32]);
                            ymm1 = _mm256_loadu_pd(&xim[32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[72],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[40]);
                            ymm1 = _mm256_loadu_pd(&xim[40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[88],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[48]);
                            ymm1 = _mm256_loadu_pd(&xim[48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[104],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[56]);
                            ymm1 = _mm256_loadu_pd(&xim[56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[120],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[64]);
                            ymm1 = _mm256_loadu_pd(&xim[64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[136],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[72]);
                            ymm1 = _mm256_loadu_pd(&xim[72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[152],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[80]);
                            ymm1 = _mm256_loadu_pd(&xim[80]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[168],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[88]);
                            ymm1 = _mm256_loadu_pd(&xim[88]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[184],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[96]);
                            ymm1 = _mm256_loadu_pd(&xim[96]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[200],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[104]);
                            ymm1 = _mm256_loadu_pd(&xim[104]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[216],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[112]);
                            ymm1 = _mm256_loadu_pd(&xim[112]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[232],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[120]);
                            ymm1 = _mm256_loadu_pd(&xim[120]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[248],ymm5);
                       }
                      else if(n > 128) {
                            register __m256d ymm0,ymm1,ymm2,ymm3;
                            register __m256d ymm4,ymm5;
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
                            ymm0 = _mm256_loadu_pd(&xre[i+0]);
                            ymm1 = _mm256_loadu_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+8]);
                            ymm1 = _mm256_loadu_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+16]);
                            ymm1 = _mm256_loadu_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+24]);
                            ymm1 = _mm256_loadu_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+56],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+32]);
                            ymm1 = _mm256_loadu_pd(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+72],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+40]);
                            ymm1 = _mm256_loadu_pd(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+88],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+48]);
                            ymm1 = _mm256_loadu_pd(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+104],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+56]);
                            ymm1 = _mm256_loadu_pd(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+120],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+64]);
                            ymm1 = _mm256_loadu_pd(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+136],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+72]);
                            ymm1 = _mm256_loadu_pd(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+152],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+80]);
                            ymm1 = _mm256_loadu_pd(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+168],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+88]);
                            ymm1 = _mm256_loadu_pd(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+184],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+96]);
                            ymm1 = _mm256_loadu_pd(&xim[i+96]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+192],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+200],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+104]);
                            ymm1 = _mm256_loadu_pd(&xim[i+104]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+208],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+216],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+112]);
                            ymm1 = _mm256_loadu_pd(&xim[i+112]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+224],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+232],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+120]);
                            ymm1 = _mm256_loadu_pd(&xim[i+120]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+240],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+248],ymm5);
                     }

                       for(; (i+95) < n; i += 96) {
                            ymm0 = _mm256_loadu_pd(&xre[i+0]);
                            ymm1 = _mm256_loadu_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+8]);
                            ymm1 = _mm256_loadu_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+16]);
                            ymm1 = _mm256_loadu_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+24]);
                            ymm1 = _mm256_loadu_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+56],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+32]);
                            ymm1 = _mm256_loadu_pd(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+72],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+40]);
                            ymm1 = _mm256_loadu_pd(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+88],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+48]);
                            ymm1 = _mm256_loadu_pd(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+104],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+56]);
                            ymm1 = _mm256_loadu_pd(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+120],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+64]);
                            ymm1 = _mm256_loadu_pd(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_store_pd(&pyc[i+136],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+72]);
                            ymm1 = _mm256_loadu_pd(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+152],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+80]);
                            ymm1 = _mm256_loadu_pd(&xim[i+80]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+160],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+168],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+88]);
                            ymm1 = _mm256_loadu_pd(&xim[i+88]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+176],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+184],ymm5);
                     }    

                       for(; (i+79) < n; i += 80) {
                            ymm0 = _mm256_loadu_pd(&xre[i+0]);
                            ymm1 = _mm256_loadu_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+8]);
                            ymm1 = _mm256_loadu_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+16]);
                            ymm1 = _mm256_loadu_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+24]);
                            ymm1 = _mm256_loadu_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+56],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+32]);
                            ymm1 = _mm256_loadu_pd(&xim[i+32]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+64],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+72],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+40]);
                            ymm1 = _mm256_loadu_pd(&xim[i+40]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+80],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+88],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+48]);
                            ymm1 = _mm256_loadu_pd(&xim[i+48]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+96],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+104],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+56]);
                            ymm1 = _mm256_loadu_pd(&xim[i+56]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+112],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+120],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+64]);
                            ymm1 = _mm256_loadu_pd(&xim[i+64]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+128],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+136],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+72]);
                            ymm1 = _mm256_loadu_pd(&xim[i+72]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+144],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+152],ymm5);
                     }

                       for(; (i+31) < n; i += 32) {
                            ymm0 = _mm256_loadu_pd(&xre[i+0]);
                            ymm1 = _mm256_loadu_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+8],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+8]);
                            ymm1 = _mm256_loadu_pd(&xim[i+8]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+16],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+24],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+16]);
                            ymm1 = _mm256_loadu_pd(&xim[i+16]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+32],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+40],ymm5);
                            ymm0 = _mm256_loadu_pd(&xre[i+24]);
                            ymm1 = _mm256_loadu_pd(&xim[i+24]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_storeu_pd(&pyc[i+48],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+56],ymm5);
                     }

                     for(; (i+15) < n; i += 16) {
                            ymm0 = _mm256_loadu_pd(&xre[i+0]);
                            ymm1 = _mm256_loadu_pd(&xim[i+0]);
                            ymm2 = _mm256_unpacklo_pd(ymm0,ymm1);
                            ymm3 = _mm256_unpackhi_pd(ymm0,ymm1);
                            ymm4 = _mm256_shuffle_f32x4(ymm2,ymm3,0x0);
                            _mm256_store_pd(&pyc[i+0],ymm4);
                            ymm5 = _mm256_shuffle_f32x4(ymm2,ymm3,0x3);
                            _mm256_storeu_pd(&pyc[i+8],ymm5);
                           
                     }

                       for(; (i+0) < n; i += 1) {
                             const double re = xre[i];
                             const double im = xim[i];
                             pyc[i]         = re;
                             pyc[i+1]       = im;
                      
                       }

                }


 
            }
   



      } // math


} // gms















#endif /*__GMS_COMPLEX_YMM4R8_HPP__*/
