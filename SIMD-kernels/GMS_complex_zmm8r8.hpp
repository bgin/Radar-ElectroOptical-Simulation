
#ifndef __GMS_COMPLEX_ZMM8R8_HPP__
#define __GMS_COMPLEX_ZMM8R8_HPP__ 261220221702


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

    const unsigned int GMS_COMPLEX_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_COMPLEX_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_COMPLEX_ZMM8R8_FULLVER =
      1000U*GMS_COMPLEX_ZMM8R8_MAJOR+
      100U*GMS_COMPLEX_ZMM8R8_MINOR+
      10U*GMS_COMPLEX_ZMM8R8_MICRO;
    const char * const GMS_COMPLEX_ZMM8R8_CREATION_DATE = "26-12-2022 17:02  +00200 (MON 26 DEC 2022 GMT+2)";
    const char * const GMS_COMPLEX_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_ZMM8R8_DESCRIPTION   = "AVX512 optimized complex number implementation."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace  gms {


       namespace math {
       
       
                      struct __ATTR_ALIGN__(64) zmm8c8_t {
                   
                          __m512d re;
                          __m512d im;
                   };

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm8r8_u( const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m512ddddd zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_pd(&xre[0]);
                        zmm1  = _mm512_loadu_pd(&yre[0]);
                        _mm512_storeu_pd(&zre[0], _mm512_add_pd(zmm0,zmm1));
                        zmm2  = _mm512_loadu_pd(&xim[0]);
                        zmm3  = _mm512_loadu_pd(&yim[0]);
                        _mm512_storeu_pd(&zim[0], _mm512_add_pd(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       const double * __restrict __ATTR_ALIGN__(64) yre,
                                       const double * __restrict __ATTR_ALIGN__(64) yim,
                                       double *       __restrict __ATTR_ALIGN__(64) zre,
                                       double *       __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512ddddd zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_pd(&xre[0]);
                        zmm1  = _mm512_load_pd(&yre[0]);
                        _mm512_store_pd(&zre[0], _mm512_add_pd(zmm0,zmm1));
                        zmm2  = _mm512_load_pd(&xim[0]);
                        zmm3  = _mm512_load_pd(&yim[0]);
                        _mm512_store_pd(&zim[0], _mm512_add_pd(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d yre,
                                     const __m512d yim,
                                     __m512d &     zre,
                                     __m512d &     zim) {
                     
                        register __m512d zmm0,zmm1;
                        zmm0 = _mm512_add_pd(xre,yre);
                        zre  = zmm0;
                        zmm1 = _mm512_add_pd(xim,yim);
                        zim  = zmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cadd_zmm8r8( const __m512d xre,
                                         const __m512d xim,
                                         const __m512d yre,
                                         const __m512d yim) {
                                     
                        zmm8c8_t cv;
                        cv.re = _mm512_add_pd(xre,yre);
                        cv.im = _mm512_add_pd(xim,yim);
                        return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cadd_zmm8r8( const zmm8c8_t x,
                                         const zmm8c8_t y) {
                                     
                        zmm8c8_t cv;
                        cv.re = _mm512_add_pd(x.re,y.re);
                        cv.im = _mm512_add_pd(x.im,y.im);
                        return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d s,
                                     __m512d &     zre,
                                     __m512d &     zim) {

                        zre = _mm512_add_pd(xre,s);
                        zim = _mm512_add_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cadd_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d s) {
                                 
                        zmm8c8_t cv;
                        cv.re = _mm512_add_pd(xre,s);
                        cv.im = _mm512_add_pd(xim,s);
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cadd_zmm8r8(const zmm8c8_t x,
                                     const __m512d s) {
                                 
                        zmm8c8_t cv;
                        cv.re = _mm512_add_pd(x.re,s);
                        cv.im = _mm512_add_pd(x.im,s);
                        return (cv);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm8r8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {
                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_pd(&xre[0]);
                        zmm1  = _mm512_loadu_pd(&xim[0]);
                        zmm2  = _mm512_loadu_pd(&zre[0]);
                        zmm3  = _mm512_loadu_pd(&zim[0])
                        _mm512_storeu_pd(&zre[0], _mm512_add_pd(zmm2,zmm0));
                        _mm512_storeu_pd(&zim[0], _mm512_add_pd(zmm3,zmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm8r8_aip(const double * __restrict __ATTR_ALIGN__(64) xre,
                                         const double * __restrict __ATTR_ALIGN__(64) xim,
                                         double *       __restrict __ATTR_ALIGN__(64) zre,
                                         double *       __restrict __ATTR_ALIGN__(64) zim) {
                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_pd(&xre[0]);
                        zmm1  = _mm512_load_pd(&xim[0]);
                        zmm2  = _mm512_load_pd(&zre[0]);
                        zmm3  = _mm512_load_pd(&zim[0])
                        _mm512_store_pd(&zre[0], _mm512_add_pd(zmm2,zmm0));
                        _mm512_store_pd(&zim[0], _mm512_add_pd(zmm3,zmm1));
              }


                ////////////////////////////////////////////////////////////////////

                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm8r8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_pd(&xre[0]);
                        zmm1  = _mm512_loadu_pd(&yre[0]);
                        _mm512_storeu_pd(&zre[0], _mm512_sub_pd(zmm0,zmm1));
                        zmm2  = _mm512_loadu_pd(&xim[0]);
                        zmm3  = _mm512_loadu_pd(&yim[0]);
                        _mm512_storeu_pd(&zim[0], _mm512_sub_pd(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       const double * __restrict __ATTR_ALIGN__(64) yre,
                                       const double * __restrict __ATTR_ALIGN__(64) yim,
                                       double *       __restrict __ATTR_ALIGN__(64) zre,
                                       double *       __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_pd(&xre[0]);
                        zmm1  = _mm512_load_pd(&yre[0]);
                        _mm512_store_pd(&zre[0], _mm512_sub_pd(zmm0,zmm1));
                        zmm2  = _mm512_load_pd(&xim[0]);
                        zmm3  = _mm512_load_pd(&yim[0]);
                        _mm512_store_pd(&zim[0], _mm512_sub_pd(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d yre,
                                     const __m512d yim,
                                     __m512d &     zre,
                                     __m512d &     zim) {
                     
                        register __m512d zmm0,zmm1;
                        zmm0 = _mm512_sub_pd(xre,yre);
                        zre  = zmm0;
                        zmm1 = _mm512_sub_pd(xim,yim);
                        zim  = zmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t csub_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d yre,
                                     const __m512d yim) {
                                    
                        zmm8c8_t cv;
                        cv.re = _mm512_sub_pd(xre,yre);
                        cv.im = _mm512_sub_pd(xim,yim);
                        return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t csub_zmm8r8(const zmm8c8_t x,
                                        const zmm8c8_t y) {
                                    
                        zmm8c8_t cv;
                        cv.re = _mm512_sub_pd(x.re,y.re);
                        cv.im = _mm512_sub_pd(x.im,y.im);
                        return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d s,
                                     __m512d &     zre,
                                     __m512d &     zim) {

                        zre = _mm512_sub_pd(xre,s);
                        zim = _mm512_sub_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t csub_zmm8r8(const __m512d xre,
                                        const __m512d xim,
                                        const __m512d s) {
                                    
                        zmm8c8_t cv;
                        cv.re = _mm512_sub_pd(xre,s);
                        cv.im = _mm512_sub_pd(xim,s);
                        return (cv);
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t csub_zmm8r8(const zmm8c8_t x,
                                        const __m512d s) {
                                    
                        zmm8c8_t cv;
                        cv.re = _mm512_sub_pd(x.re,s);
                        cv.im = _mm512_sub_pd(x.im,s);
                        return (cv);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm8r8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {
                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_pd(&xre[0]);
                        zmm1  = _mm512_loadu_pd(&xim[0]);
                        zmm2  = _mm512_loadu_pd(&zre[0]);
                        zmm3  = _mm512_loadu_pd(&zim[0])
                        _mm512_storeu_pd(&zre[0], _mm512_sub_pd(zmm2,zmm0));
                        _mm512_storeu_pd(&zim[0], _mm512_sub_pd(zmm3,zmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm8r8_aip(const double * __restrict __ATTR_ALIGN__(64) xre,
                                         const double * __restrict __ATTR_ALIGN__(64) xim,
                                         double *       __restrict __ATTR_ALIGN__(64) zre,
                                         double *       __restrict __ATTR_ALIGN__(64) zim) {
                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_pd(&xre[0]);
                        zmm1  = _mm512_load_pd(&xim[0]);
                        zmm2  = _mm512_load_pd(&zre[0]);
                        zmm3  = _mm512_load_pd(&zim[0])
                        _mm512_store_pd(&zre[0], _mm512_sub_pd(zmm2,zmm0));
                        _mm512_store_pd(&zim[0], _mm512_sub_pd(zmm3,zmm1));
              }


               ////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm8r8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                           register __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_loadu_pd(&xre[0]);
                           zmm1  = _mm512_loadu_pd(&yre[0]);
                           zmm2  = _mm512_loadu_pd(&xim[0]);
                           zmm3  = _mm512_loadu_pd(&yim[0]);
                           zmm4  = _mm512_sub_pd(_mm512_mul_pd(zmm0,zmm1),
                                                                        _mm512_mul_pd(zmm2,zmm3));
                           _mm512_storeu_pd(&zre[0], zmm4);
                           zmm5  = _mm512_mul_pd(_mm512_mul_pd(zmm2,zmm1),
                                                                        _mm512_mul_pd(zmm0,zmm3));
                           _mm512_storeu_pd(&zim[0], zmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       const double * __restrict __ATTR_ALIGN__(64) yre,
                                       const double * __restrict __ATTR_ALIGN__(64) yim,
                                       double *       __restrict __ATTR_ALIGN__(64) zre,
                                       double *       __restrict __ATTR_ALIGN__(64) zim) {

                           register __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_load_pd(&xre[0]);
                           zmm1  = _mm512_load_pd(&yre[0]);
                           zmm2  = _mm512_load_pd(&xim[0]);
                           zmm3  = _mm512_load_pd(&yim[0]);
                           zmm4  = _mm512_sub_pd(_mm512_mul_pd(zmm0,zmm1),
                                                                        _mm512_mul_pd(zmm2,zmm3));
                           _mm512_store_pd(&zre[0], zmm4);
                           zmm5  = _mm512_mul_pd(_mm512_mul_pd(zmm2,zmm1),
                                                                        _mm512_mul_pd(zmm0,zmm3));
                           _mm512_store_pd(&zim[0], zmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d yre,
                                     const __m512d yim,
                                     __m512d &     zre,
                                     __m512d &     zim) {

                         register __m512d zmm0,zmm1;
                         zmm0 = _mm512_sub_pd(_mm512_mul_pd(xre,yre),
                                              _mm512_mul_pd(xim,yim));
                         zre  = zmm0;
                         zmm1 = _mm512_mul_pd(_mm512_mul_pd(xim,yre),
                                              _mm512_mul_pd(xre,yim));
                         zim  = zmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cmul_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d yre,
                                     const __m512d yim) {
                                     

                         zmm8c8_t cv;
                         cv.re = _mm512_sub_pd(_mm512_mul_pd(xre,yre),
                                              _mm512_mul_pd(xim,yim));
                         cv.im = _mm512_mul_pd(_mm512_mul_pd(xim,yre),
                                              _mm512_mul_pd(xre,yim));
                         return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cmul_zmm8r8(const zmm8c8_t x,
                                        const zmm8c8_t y) {
                                     

                         zmm8c8_t cv;
                         cv.re = _mm512_sub_pd(_mm512_mul_pd(x.re,y.re),
                                              _mm512_mul_pd(x.im,y.im));
                         cv.im = _mm512_mul_pd(_mm512_mul_pd(x.im,y.re),
                                              _mm512_mul_pd(x.re,y.im));
                         return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     const __m512d s,
                                     __m512d &     zre,
                                     __m512d &     zim) {

                        zre = _mm512_mul_pd(xre,s);
                        zim = _mm512_mul_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cmul_zmm8r8( const zmm8c8_t x
                                     const __m512d s) {
                                    
                        zmm8c8_t cv;
                        zre = _mm512_mul_pd(x.re,s);
                        zim = _mm512_mul_pd(x.im,s);
                        return (cv);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm8r8_uip(const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {

                           register __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_loadu_pd(&xre[0]);
                           zmm1  = _mm512_loadu_pd(&zre[0]);
                           zmm2  = _mm512_loadu_pd(&xim[0]);
                           zmm3  = _mm512_loadu_pd(&zim[0]);
                           zmm4  = _mm512_sub_pd(_mm512_mul_pd(zmm0,zmm1),
                                                 _mm512_mul_pd(zmm2,zmm3));
                           _mm512_storeu_pd(&zre[0], zmm4);
                           zmm5  = _mm512_mul_pd(_mm512_mul_pd(zmm2,zmm1),
                                                 _mm512_mul_pd(zmm0,zmm3));
                           _mm512_storeu_pd(&zim[0], zmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm8r8_aip(const double * __restrict __ATTR_ALIGN__(64) xre,
                                         const double * __restrict __ATTR_ALIGN__(64) xim,
                                         double *       __restrict __ATTR_ALIGN__(64) zre,
                                         double *       __restrict __ATTR_ALIGN__(64) zim) {

                           register __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_load_pd(&xre[0]);
                           zmm1  = _mm512_load_pd(&zre[0]);
                           zmm2  = _mm512_load_pd(&xim[0]);
                           zmm3  = _mm512_load_pd(&zim[0]);
                           zmm4  = _mm512_sub_pd(_mm512_mul_pd(zmm0,zmm1),
                                                 _mm512_mul_pd(zmm2,zmm3));
                           _mm512_store_pd(&zre[0], zmm4);
                           zmm5  = _mm512_mul_pd(_mm512_mul_pd(zmm2,zmm1),
                                                 _mm512_mul_pd(zmm0,zmm3));
                           _mm512_store_pd(&zim[0], zmm5);
               }

                 ////////////////////////////////////////////////////////////////////////////


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8_u( const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3; 
                        register __m512d zmm4,zmm5,zmm6;
                        zmm0  = _mm512_loadu_pd(&xre[0]); //a
                        zmm1  = _mm512_loadu_pd(&yim[0]); //d
                        zmm2  = _mm512_loadu_pd(&xim[0]); //b
                        zmm3  = _mm512_loadu_pd(&yre[0]); //c
                        zmm4  = _mm512_fmadd_pd(zmm0,zmm3,
                                                _mm512_mul_pd(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_pd(zmm2,zmm3,
                                                _mm512_mul_pd(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_pd(zmm3,zmm3),
                                                _mm512_mul_pd(zmm1,zmm1));
                        _mm512_storeu_pd(&zre[0], _mm512_div_pd(zmm4,zmm6));
                        _mm512_storeu_pd(&zim[0], _mm512_div_pd(zmm5,zmm6));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8_a( const double * __restrict xre,
                                       const double * __restrict xim,
                                       const double * __restrict yre,
                                       const double * __restrict yim,
                                       double *       __restrict zre,
                                       double *       __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3; 
                        register __m512d zmm4,zmm5,zmm6;
                        zmm0  = _mm512_load_pd(&xre[0]); //a
                        zmm1  = _mm512_load_pd(&yim[0]); //d
                        zmm2  = _mm512_load_pd(&xim[0]); //b
                        zmm3  = _mm512_load_pd(&yre[0]); //c
                        zmm4  = _mm512_fmadd_pd(zmm0,zmm3,
                                                _mm512_mul_pd(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_pd(zmm2,zmm3,
                                                _mm512_mul_pd(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_pd(zmm3,zmm3,
                                                _mm512_mul_pd(zmm1,zmm1));
                        _mm512_store_pd(&zre[0], _mm512_div_pd(zmm4,zmm6));
                        _mm512_store_pd(&zim[0], _mm512_div_pd(zmm5,zmm6));
                }
              
                                         
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8( const __m512d xre,
                                     const __m512d xim,
                                     const __m512d yre,
                                     const __m512d yim,
                                     __m512d & zre,
                                     __m512d & zim) {

                      register __m512d zmm0,zmm1,zmm2;
                      zmm0 = _mm512_fmadd_pd(xre,yre,
                                           _mm512_mul_pd(xim,yim));
                      zmm1 = _mm512_fmsub_pd(xim,yre,
                                           _mm512_mul_pd(xre,yim));
                      zmm2 = _mm512_fmadd_pd(zmm3,zmm3,
                                           _mm512_mul_pd(zmm1,zmm1));
                      zre  = _mm512_div_pd(zmm0,zmm2);
                      zim  = _mm512_div_pd(zmm1,zmm2);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_zmm8r8( const __m512d xre,
                                         const __m512d xim,
                                         const __m512d yre,
                                         const __m512d yim) {
                         
                      zmm8c8_t vc;           
                      register __m512d zmm0,zmm1,zmm2;
                      zmm0 = _mm512_fmadd_pd(xre,yre,
                                           _mm512_mul_pd(xim,yim));
                      zmm1 = _mm512_fmsub_pd(xim,yre,
                                           _mm512_mul_pd(xre,yim));
                      zmm2 = _mm512_fmadd_pd(zmm3,zmm3,
                                           _mm512_mul_pd(zmm1,zmm1));
                      cv.re  = _mm512_div_pd(zmm0,zmm2);
                      cv.im  = _mm512_div_pd(zmm1,zmm2);
                      return (cv);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_zmm8r8( const zmm8c8_t x,
                                         const zmm8c8_t y) {
                         
                      zmm8c8_t vc;           
                      register __m512d zmm0,zmm1,zmm2;
                      zmm0 = _mm512_fmadd_pd(x.re,y.re,
                                           _mm512_mul_pd(x.im,y.im));
                      zmm1 = _mm512_fmsub_pd(x.im,y.re,
                                           _mm512_mul_pd(x.re,y.im));
                      zmm2 = _mm512_fmadd_pd(zmm3,zmm3,
                                           _mm512_mul_pd(zmm1,zmm1));
                      cv.re  = _mm512_div_pd(zmm0,zmm2);
                      cv.im  = _mm512_div_pd(zmm1,zmm2);
                      return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8( const __m512d xre,
                                     const __m512d xim,
                                     const __m512d s,
                                     __m512d & zre,
                                     __m512d & zim) {

                        zre = _mm512_div_pd(xre,s);
                        zim = _mm512_div_pd(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_zmm8r8( const __m512d xre,
                                     const __m512d xim,
                                     const __m512d s) {
                                     
                        zmm8c8_t cv;
                        cv.re = _mm512_div_pd(xre,s);
                        cv.im = _mm512_div_pd(xim,s);
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_zmm8r8( const zmm8c8_t x,
                                     const __m512d s) {
                                     
                        zmm8c8_t cv;
                        cv.re = _mm512_div_pd(x.re,s);
                        cv.im = _mm512_div_pd(x.im,s);
                        return (cv);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8_uip( const double * __restrict xre,
                                         const double * __restrict xim,
                                         double *       __restrict zre,
                                         double *       __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3; 
                        register __m512d zmm4,zmm5,zmm6;
                        zmm0  = _mm512_loadu_pd(&xre[0]); //a
                        zmm1  = _mm512_loadu_pd(&zim[0]); //d
                        zmm2  = _mm512_loadu_pd(&xim[0]); //b
                        zmm3  = _mm512_loadu_pd(&zre[0]); //c
                        zmm4  = _mm512_fmadd_pd(zmm0,zmm3,
                                                _mm512_mul_pd(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_pd(zmm2,zmm3,
                                                _mm512_mul_pd(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_pd(zmm3,zmm3,
                                                _mm512_mul_pd(zmm1,zmm1));
                        _mm512_storeu_pd(&zre[0], _mm512_div_pd(zmm4,zmm6));
                        _mm512_storeu_pd(&zim[0], _mm512_div_pd(zmm5,zmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8_aip( const double * __restrict __ATTR_ALIGN__(64) xre,
                                         const double * __restrict __ATTR_ALIGN__(64) xim,
                                         double *       __restrict __ATTR_ALIGN__(64) zre,
                                         double *       __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3; 
                        register __m512d zmm4,zmm5,zmm6;
                        zmm0  = _mm512_load_pd(&xre[0]); //a
                        zmm1  = _mm512_load_pd(&zim[0]); //d
                        zmm2  = _mm512_load_pd(&xim[0]); //b
                        zmm3  = _mm512_load_pd(&zre[0]); //c
                        zmm4  = _mm512_fmadd_pd(zmm0,zmm3,
                                                _mm512_mul_pd(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_pd(zmm2,zmm3,
                                                _mm512_mul_pd(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_pd(zmm3,zmm3,
                                                _mm512_mul_pd(zmm1,zmm1));
                        _mm512_store_pd(&zre[0], _mm512_div_pd(zmm4,zmm6));
                        _mm512_store_pd(&zim[0], _mm512_div_pd(zmm5,zmm6));
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm8r8_s( const __m512d s,
                                       const __m512d xre,
                                       const __m512d xim,
                                       __m512d * __restrict zre,
                                       __m512d * __restrict zim) {
                        
                        register __m512d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_pd();
                        cdiv_zmm8r8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_zmm8r8_s(const __m512d s,
                                            const __m512d xre,
                                            const __m512d xim) {
                                       
                        zmm8c8_t cv;
                        register __m512d t0r,t0i;
                        t0r = s;
                        t0i = _mm512_setzero_pd();
                        cdiv_zmm8r8(t0r,t0i,xre,xim,&cv.re,&cv.im);
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_zmm8r8_s(const __m512d s,
                                          const zmm8c8_t x) {
                                       
                        zmm8c8_t cv;
                        register __m512d t0r,t0i;
                        t0r = s;
                        t0i = _mm512_setzero_pd();
                        cdiv_zmm8r8(t0r,t0i,x.re,x.im,&cv.re,&cv.im);
                        return (cv);                     
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm8r8_u( const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d r,den;
                        __mmask8 m = 0x0;
                        zmm0 = _mm512_loadu_pd(&yre[0]); // c
                        zmm1 = _mm512_loadu_pd(&yim[0]); // d
                        zmm2 = _mm512_loadu_pd(&xre[0]); // a
                        zmm3 = _mm512_loadu_pd(&xim[0]); // b
                        m    = _mm512_cmp_pd_mask(_mm512_abs_pd(zmm0),
                                                  _mm512_abs_pd(zmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_pd(m,_mm512_div_pd(zmm0,zmm1),
                                                      _mm512_div_pd(zmm1,zmm0)); // r
                        den  = _mm512_mask_blend_pd(m,_mm512_fmadd_pd(r,zmm0,zmm1),
                                                      _mm512_fmadd_pd(r,zmm1,zmm0));
                        _mm512_storeu_pd(&zre[0], _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmadd_pd(zmm2,r,zmm3),den),
                                                _mm512_div_pd(_mm512_fmadd_pd(zmm3,r,zmm2),den)));
                        _mm512_storeu_pd(&zim[0], _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmsub_pd(zmm3,r,zmm2),den),
                                                _mm512_div_pd(_mm512_sub_pd(zmm3,_mm512_mul_pd(zmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm8r8_a( const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d r,den;
                        __mmask8 m = 0x0;
                        zmm0 = _mm512_load_pd(&yre[0]); // c
                        zmm1 = _mm512_load_pd(&yim[0]); // d
                        zmm2 = _mm512_load_pd(&xre[0]); // a
                        zmm3 = _mm512_load_pd(&xim[0]); // b
                        m    = _mm512_cmp_pd_mask(_mm512_abs_pd(zmm0),
                                                  _mm512_abs_pd(zmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_pd(m,_mm512_div_pd(zmm0,zmm1),
                                                      _mm512_div_pd(zmm1,zmm0)); // r
                        den  = _mm512_mask_blend_pd(m,_mm512_fmadd_pd(r,zmm0,zmm1),
                                                      _mm512_fmadd_pd(r,zmm1,zmm0));
                        _mm512_storeu_pd(&zre[0], _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmadd_pd(zmm2,r,zmm3),den),
                                                _mm512_div_pd(_mm512_fmadd_pd(zmm3,r,zmm2),den)));
                        _mm512_storeu_pd(&zim[0], _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmsub_pd(zmm3,r,zmm2),den),
                                                _mm512_div_pd(_mm512_sub_pd(zmm3,_mm512_mul_pd(zmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm8r8( const __m512d xre,
                                           const __m512d xim,
                                           const __m512d yre,
                                           const __m512d yim,
                                           __m512d * __restrict zre,
                                           __m512d * __restrict zim) {

                        register __m512d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm512_cmp_pd_mask(_mm512_abs_pd(yre),
                                                  _mm512_abs_pd(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_pd(m,_mm512_div_pd(yre,yim),
                                                      _mm512_div_pd(yim,yre)); // r
                        den  = _mm512_mask_blend_pd(m,_mm512_fmadd_pd(r,yre,yim),
                                                      _mm512_fmadd_pd(r,yim,yre));
                        *zre  =  _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmadd_pd(xre,r,xim),den),
                                                _mm512_div_pd(_mm512_fmadd_pd(xim,r,xre),den));
                        *zim  =  _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmsub_pd(xim,r,xre),den),
                                                _mm512_div_pd(_mm512_sub_pd(xim,_mm512_mul_pd(xre,r)),den)));
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_smith_zmm8r8( const __m512d xre,
                                           const __m512d xim,
                                           const __m512d yre,
                                           const __m512d yim) {
                                          

                        zmm8c8_t cv;
                        register __m512d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm512_cmp_pd_mask(_mm512_abs_pd(yre),
                                                  _mm512_abs_pd(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_pd(m,_mm512_div_pd(yre,yim),
                                                      _mm512_div_pd(yim,yre)); // r
                        den  = _mm512_mask_blend_pd(m,_mm512_fmadd_pd(r,yre,yim),
                                                      _mm512_fmadd_pd(r,yim,yre));
                        cv.re  =  _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmadd_pd(xre,r,xim),den),
                                                _mm512_div_pd(_mm512_fmadd_pd(xim,r,xre),den));
                        cv.im  =  _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmsub_pd(xim,r,xre),den),
                                                _mm512_div_pd(_mm512_sub_pd(xim,_mm512_mul_pd(xre,r)),den)));
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_smith_zmm8r8( const zmm8c8_t x,
                                               const zmm8c8_t y) {
                                          

                        zmm8c8_t cv;
                        register __m512d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm512_cmp_pd_mask(_mm512_abs_pd(y.re),
                                                  _mm512_abs_pd(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_pd(m,_mm512_div_pd(y.re,y.im),
                                                      _mm512_div_pd(y.im,y.re)); // r
                        den  = _mm512_mask_blend_pd(m,_mm512_fmadd_pd(r,y.re,yim),
                                                      _mm512_fmadd_pd(r,y.im,y.re));
                        cv.re  =  _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmadd_pd(x.re,r,x.im),den),
                                                _mm512_div_pd(_mm512_fmadd_pd(x.im,r,x.re),den));
                        cv.im  =  _mm512_mask_blend_pd(m,
                                                _mm512_div_pd(_mm512_fmsub_pd(x.im,r,x.re),den),
                                                _mm512_div_pd(_mm512_sub_pd(x.im,_mm512_mul_pd(x.re,r)),den)));
                        return (cv);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm8r8_s(const __m512d s,
                                             const __m512d xre,
                                             const __m512d xim,
                                             __m512d * __restrict zre,
                                             __m512d * __restrict zim) {
                                             
                        register __m512d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_pd(); 
                        cdiv_smith_zmm8r8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                   
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_smith_zmm8r8_s(const __m512d s,
                                                  const __m512d xre,
                                                  const __m512d xim) {
                                             
                        zmm8c8_t cv;                    
                        register __m512d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_pd(); 
                        cdiv_smith_zmm8r8(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm8c8_t cdiv_smith_zmm8r8_s(const __m512d s,
                                                const zmm8c8_t x) {
                                             
                        zmm8c8_t cv;                    
                        register __m512d t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_pd(); 
                        cdiv_smith_zmm8r8(t0r,t0i,x.re,x.im,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
               


#include "GMS_sleefsimddp.hpp"


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_zmm8r8_u(const double * __restrict re,
                                       const double * __restrict im,
                                       double * __restrict  cabs) {

                        register __m512d zmm0,zmm1,zmm2,zmm3,zmm4;
                        zmm0  = _mm512_loadu_pd(&re[0]);
                        zmm1  = _mm512_mul_pd(zmm0,zmm0);
                        zmm2  = _mm512_loadu_pd(&im[0]);
                        zmm3  = _mm512_mul_pd(zmm2,zmm2);
                        zmm4  = xsqrt(_mm512_add_pd(zmm1,zmm3));
                        _mm512_storeu_pd(&cabs[0],zmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) re,
                                       const double * __restrict __ATTR_ALIGN__(64) im,
                                       double * __restrict  __ATTR_ALIGN__(64) cabs) {

                        register __m512d zmm0,zmm1,zmm2,zmm3,zmm4;
                        zmm0  = _mm512_load_pd(&re[0]);
                        zmm1  = _mm512_mul_pd(zmm0,zmm0);
                        zmm2  = _mm512_load_pd(&im[0]);
                        zmm3  = _mm512_mul_pd(zmm2,zmm2);
                        zmm4  = xsqrt(_mm512_add_pd(zmm1,zmm3));
                        _mm512_store_pd(&cabs[0],zmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d cabs_zmm8r8(const __m512d re,
                                       const __m512d im) {

                        register __m512d zmm0,zmm1,cabs;
                        zmm0 = _mm512_mul_pd(re,re);
                        zmm1 = _mm512_mul_pd(im,im);
                        cabs = xsqrt(_mm512_add_pd(zmm0,zmm1));
                        return (cabs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d cabs_zmm8r8(const zmm8c8_t x) {

                        register __m512d zmm0,zmm1,cabs;
                        zmm0 = _mm512_mul_pd(x.re,x.re);
                        zmm1 = _mm512_mul_pd(x.im,x.im);
                        cabs = xsqrt(_mm512_add_pd(zmm0,zmm1));
                        return (cabs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_zmm8r8_u(const double * __restrict re,
                                       const double * __restrict im,
                                       double * __restrict  carg) {

                        register __m512d zmm0,zmm1;
                        zmm0 = _mm512_loadu_pd(&re[0]);
                        zmm1 = _mm512_loadu_pd(&im[0]);
                        _mm512_storeu_pd(&carg[0], xatan2(zmm0,zmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) re,
                                       const double * __restrict __ATTR_ALIGN__(64) im,
                                       double * __restrict  __ATTR_ALIGN__(64) carg) {

                        register __m512d zmm0,zmm1;
                        zmm0 = _mm512_load_pd(&re[0]);
                        zmm1 = _mm512_load_pd(&im[0]);
                        _mm512_store_pd(&carg[0], xatan2(zmm0,zmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d carg_zmm8r8(const __m512d re,
                                       const __m512d im) {

                       register __m512d carg;
                       carg = xatan2(re,im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d carg_zmm8r8(const zmm8c8_t x) {

                       register __m512d carg;
                       carg = xatan2(x.re,x.im);
                       return (carg);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm8r8_v2(const __m512d xre,
                                         const __m512d xim,
                                         __m512d * __restrict yre,
                                         __m512d * __restrict yim) {
                         
                        //register __m512 c;              
                        //c = negate_zmm16r4(*im);
                        //*im = c;
                        *yre = xre; 
                        *yim = negate_zmm8r8(xim);
                   } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm8r8_u( double * __restrict re,
                                        double * __restrict  im) {

                        register __m512d c;
                        c = negate_zmm8r8(_mm512_loadu_pd(&im[0]));
                        _mm512_storeu_pd(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64)  re,
                                       double * __restrict  __ATTR_ALIGN__(64) im) {

                        register __m512d c;
                        c = negate_zmm8r8(_mm512_load_pd(&im[0]));
                        _mm512_store_pd(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512ddd cconj_zmm8r8( __m512d * __restrict re,
                                           __m512d * __restrict im ) {
                          
                        register __m512d c;
                        c = negate_zmm8r8(*im);
                        *im = c;
               } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_zmm8r8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict  csre,
                                       double * __restrict  csim) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_loadu_pd(&xre[0]);
                      zmm1  = _mm512_loadu_pd(&xim[0]);
                      zmm2  = _mm512_mul_pd(xcos(zmm0),xcosh(zmm1));
                      _mm512_storeu_pd(&csre[0],zmm2);
                      zmm3  = _mm512_mul_pd(xsin(zmm0),xsinh(zmm1));
                      _mm512_storeu_pd(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       double * __restrict  __ATTR_ALIGN__(64) csre,
                                       double * __restrict  __ATTR_ALIGN__(64) csim) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_load_pd(&xre[0]);
                      zmm1  = _mm512_load_pd(&xim[0]);
                      zmm2  = _mm512_mul_pd(xcos(zmm0),xcosh(zmm1));
                      _mm512_store_pd(&csre[0],zmm2);
                      zmm3  = _mm512_mul_pd(xsin(zmm0),xsinh(zmm1));
                      _mm512_store_pd(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     __m512d * __restrict csre,
                                     __m512d * __restrict csim) {

                      register __m512d zmm0,zmm1;
                      zmm0  = _mm512_mul_pd(xcos(xre),xcosh(xim));
                      *csre = zmm0;
                      zmm1  = _mm512_mul_pd(xsin(xre),xsinh(xim));
                      *csim = zmm1; 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_zmm8r8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict  csre,
                                       double * __restrict  csim) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_loadu_pd(&xre[0]);
                      zmm1  = _mm512_loadu_pd(&xim[0]);
                      zmm2  = _mm512_mul_pd(xcosh(zmm0),xcos(zmm1));
                      _mm512_storeu_pd(&csre[0],zmm2);
                      zmm3  = _mm512_mul_pd(xsinh(zmm0),xsin(zmm1));
                      _mm512_storeu_pd(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       double * __restrict  __ATTR_ALIGN__(64) csre,
                                       double * __restrict  __ATTR_ALIGN__(64) csim) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_load_pd(&xre[0]);
                      zmm1  = _mm512_load_pd(&xim[0]);
                      zmm2  = _mm512_mul_pd(xcosh(zmm0),xcos(zmm1));
                      _mm512_store_pd(&csre[0],zmm2);
                      zmm3  = _mm512_mul_pd(xsinh(zmm0),xsin(zmm1));
                      _mm512_store_pd(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     __m512d * __restrict csre,
                                     __m512d * __restrict csim) {

                      register __m512d zmm0,zmm1;
                      zmm0  = _mm512_mul_pd(xcosh(xre),xcos(xim));
                      *csre = zmm0;
                      zmm1  = _mm512_mul_pd(xsinh(xre),xsin(xim));
                      *csim = zmm1; 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_zmm8r8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_pd(&xre[0]);
                      zmm1 = _mm512_loadu_pd(&yre[0]);
                      _mm512_storeu_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_EQ_OQ));
                      zmm2 = _mm512_loadu_pd(&xim[0]);
                      zmm3 = _mm512_loadu_pd(&yim[0]);
                      _mm512_storeu_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_EQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                      const double * __restrict __ATTR_ALIGN__(64) xim,
                                      const double * __restrict __ATTR_ALIGN__(64) yre,
                                      const double * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_pd(&xre[0]);
                      zmm1 = _mm512_load_pd(&yre[0]);
                      _mm512_store_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_EQ_OQ));
                      zmm2 = _mm512_load_pd(&xim[0]);
                      zmm3 = _mm512_load_pd(&yim[0]);
                      _mm512_store_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_EQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_zmm8r8(const __m512d xre,
                                    const __m512d xim,
                                    const __m512d yre,
                                    const __m512d yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_pd_mask(xre,yre,_CMP_EQ_OQ);
                         *eqi = _mm512_cmp_pd_mask(xim,yim,_CMP_EQ_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_zmm8r8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_pd(&xre[0]);
                      zmm1 = _mm512_loadu_pd(&yre[0]);
                      _mm512_storeu_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_GT_OQ));
                      zmm2 = _mm512_loadu_pd(&xim[0]);
                      zmm3 = _mm512_loadu_pd(&yim[0]);
                      _mm512_storeu_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_GT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                      const double * __restrict __ATTR_ALIGN__(64) xim,
                                      const double * __restrict __ATTR_ALIGN__(64) yre,
                                      const double * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_pd(&xre[0]);
                      zmm1 = _mm512_load_pd(&yre[0]);
                      _mm512_store_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_GT_OQ));
                      zmm2 = _mm512_load_pd(&xim[0]);
                      zmm3 = _mm512_load_pd(&yim[0]);
                      _mm512_store_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_GT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_zmm8r8(const __m512d xre,
                                    const __m512d xim,
                                    const __m512d yre,
                                    const __m512d yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_pd_mask(xre,yre,_CMP_GT_OQ);
                         *eqi = _mm512_cmp_pd_mask(xim,yim,_CMP_GT_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_zmm8r8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_pd(&xre[0]);
                      zmm1 = _mm512_loadu_pd(&yre[0]);
                      _mm512_storeu_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_LT_OQ));
                      zmm2 = _mm512_loadu_pd(&xim[0]);
                      zmm3 = _mm512_loadu_pd(&yim[0]);
                      _mm512_storeu_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_LT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                      const double * __restrict __ATTR_ALIGN__(64) xim,
                                      const double * __restrict __ATTR_ALIGN__(64) yre,
                                      const double * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_pd(&xre[0]);
                      zmm1 = _mm512_load_pd(&yre[0]);
                      _mm512_store_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_LT_OQ));
                      zmm2 = _mm512_load_pd(&xim[0]);
                      zmm3 = _mm512_load_pd(&yim[0]);
                      _mm512_store_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_LT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_zmm8r8(const __m512d xre,
                                    const __m512d xim,
                                    const __m512d yre,
                                    const __m512d yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_pd_mask(xre,yre,_CMP_LT_OQ);
                         *eqi = _mm512_cmp_pd_mask(xim,yim,_CMP_LT_OQ);
              }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_zmm8r8_u(const double * __restrict xre,
                                      const double * __restrict xim,
                                      const double * __restrict yre,
                                      const double * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_pd(&xre[0]);
                      zmm1 = _mm512_loadu_pd(&yre[0]);
                      _mm512_storeu_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_NEQ_OQ));
                      zmm2 = _mm512_loadu_pd(&xim[0]);
                      zmm3 = _mm512_loadu_pd(&yim[0]);
                      _mm512_storeu_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_NEQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                      const double * __restrict __ATTR_ALIGN__(64) xim,
                                      const double * __restrict __ATTR_ALIGN__(64) yre,
                                      const double * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512d zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_pd(&xre[0]);
                      zmm1 = _mm512_load_pd(&yre[0]);
                      _mm512_store_pd(&eqr[0],
                                       _mm512_cmp_pd_mask(zmm0,zmm1,_CMP_NEQ_OQ));
                      zmm2 = _mm512_load_pd(&xim[0]);
                      zmm3 = _mm512_load_pd(&yim[0]);
                      _mm512_store_pd(&eqi[0],
                                       _mm512_cmp_pd_mask(zmm2,zmm3,_CMP_NEQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_zmm8r8(const __m512d xre,
                                    const __m512d xim,
                                    const __m512d yre,
                                    const __m512d yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_pd_mask(xre,yre,_CMP_NEQ_OQ);
                         *eqi = _mm512_cmp_pd_mask(xim,yim,_CMP_NEQ_OQ);
              }


                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_zmm8r8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict cexpr,
                                       double * __restrict cexpi ) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_pd(&xre[0]);
                        zmm1  = _mm512_loadu_pd(&xim[0]);
                        zmm2  = xexp(zmm0);
                        zmm3  = _mm512_mul_pd(zmm2,xcos(zmm1));
                        _mm512_storeu_pd(&cexpr[0],zmm3);
                        zmm4  = _mm512_mul_pd(zmm2,xsin(zmm1));
                        _mm512_storeu_pd(&cexpi[0],zmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       double * __restrict __ATTR_ALIGN__(64) cexpr,
                                       double * __restrict __ATTR_ALIGN__(64) cexpi ) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_pd(&xre[0]);
                        zmm1  = _mm512_load_pd(&xim[0]);
                        zmm2  = xexp(zmm0);
                        zmm3  = _mm512_mul_pd(zmm2,xcos(zmm1));
                        _mm512_store_pd(&cexpr[0],zmm3);
                        zmm4  = _mm512_mul_pd(zmm2,xsin(zmm1));
                        _mm512_store_pd(&cexpi[0],zmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_zmm8r8(const __m512d xre,
                                     const __m512d xim,
                                     __m512d * __restrict cexpr,
                                     __m512d * __restrict cexpi) {

                        register __m512d zmm0;
                        zmm0   = xexp(xre);
                        *cexpr = _mm512_mul_pd(zmm0,xcos(xim));
                        *cexpi = _mm512_mul_pd(zmm0,xsin(xim));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_zmm8r8_u(const double * __restrict rho,
                                         const double * __restrict tht,
                                         double * __restrict  re,
                                         double * __restrict  im) {

                         register __m512d zmm0,zmm1,zmm2,zmm3;
                         zmm0 = _mm512_loadu_pd(&rho[0]);
                         zmm1 = _mm512_loadu_pd(&tht[0]);
                         zmm2 = _mm512_mul_pd(zmm0,xcos(zmm1)); //tht
                         _mm512_storeu_pd(&re[0],zmm2);
                         zmm3 = _mm512_mul_pd(zmm0,xsin(zmm1)); //tht
                         _mm512_storeu_pd(&im[0],zmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) rho,
                                         const double * __restrict __ATTR_ALIGN__(64) tht,
                                         double * __restrict  __ATTR_ALIGN__(64) re,
                                         double * __restrict  __ATTR_ALIGN__(64) im) {

                         register __m512d zmm0,zmm1,zmm2,zmm3;
                         zmm0 = _mm512_load_pd(&rho[0]);
                         zmm1 = _mm512_load_pd(&tht[0]);
                         zmm2 = _mm512_mul_pd(zmm0,xcos(zmm1)); //tht
                         _mm512_store_pd(&re[0],zmm2);
                         zmm3 = _mm512_mul_pd(zmm0,xsin(zmm1)); //tht
                         _mm512_store_pd(&im[0],zmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_zmm8r8(const __m512d rho,
                                       const __m512d tht,
                                       __m512d * __restrict re,
                                       __m512d * __restrict im) {

                        register __m512d zmm0,zmm1;
                        zmm0 = _mm512_mul_pd(rho,xcos(tht));
                        *re  = zmm0;
                        zmm1 = _mm512_mul_pd(rho,xsin(tht));
                        *im  = zmm1;
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_zmm8r8_u(const double * __restrict xre,
                                       const double * __restrict xim,
                                       double * __restrict wrkc,
                                       double * __restrict csqr,
                                       double * __restrict csqi) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        const register __m512dd half = _mm512_set1_pd(0.5);
                        cabs_zmm8r8_u(xre,xim,wrkc);
                        zmm0  = _mm512_loadu_pd(&xre[0]);
                        zmm1  = _mm512_loadu_pd(&wrkc[0]);
                        zmm2  = _mm512_mul_pd(half,_mm512_add_pd(zmm1,zmm0));
                        _mm512_storeu_pd(&csqr[0],_mm512_sqrt_pd(zmm2));
                        zmm3  = _mm512_mul_pd(half,_mm512_sub_pd(zmm1,zmm0));
                        _mm512_storeu_pd(&csqi[0],_mm512_sqrt_pd(zmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                       const double * __restrict __ATTR_ALIGN__(64) xim,
                                       double * __restrict __ATTR_ALIGN__(64) wrkc,
                                       double * __restrict __ATTR_ALIGN__(64) csqr,
                                       double * __restrict __ATTR_ALIGN__(64) csqi) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        const register __m512dd half = _mm512_set1_pd(0.5);
                        cabs_zmm8r8_a(xre,xim,wrkc);
                        zmm0  = _mm512_load_pd(&xre[0]);
                        zmm1  = _mm512_load_pd(&wrkc[0]);
                        zmm2  = _mm512_mul_pd(half,_mm512_add_pd(zmm1,zmm0));
                        _mm512_store_pd(&csqr[0],_mm512_sqrt_pd(zmm2));
                        zmm3  = _mm512_mul_pd(half,_mm512_sub_pd(zmm1,zmm0));
                        _mm512_store_pd(&csqi[0],_mm512_sqrt_pd(zmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_zmm8r8(const __m512d xre,
                                      const __m512d xim,
                                      __m512d * __restrict wrkc,
                                      __m512d * __restrict csqr,
                                      __m512d * __restrict csqi) {

                       register __m512d zmm0,zmm1;
                       const register __m512dd half = _mm512_set1_pd(0.5); 
                       cabs_zmm8r8(xre,xim,wrkc);
                       zmm0  = _mm512_mul_pd(half,_mm512_add_pd(*wrkc,xre));
                       *csqr = zmm0;
                       zmm1  = _mm512_mul_pd(half,_mm512_sub_pd(*wrkc,xre));
                       *csqi = zmm1; 
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_zmm8r8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict zre,
                                             double * __restrict zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        rep  = _mm512_fmsub_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        zmm0 = _mm512_mul_pd(rep,rep);
                        zmm1 = _mm512_mul_pd(imp,imp);
                        zmm2 = _mm512_sqrt_pd(_mm512_add_pd(zmm0,zmm1));
                        _mm512_storeu_pd(&zre[0], _mm512_div_pd(rep,zmm2));
                        _mm512_storeu_pd(&zim[0], _mm512_div_pd(imp,zmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                             const double * __restrict __ATTR_ALIGN__(64) xim,
                                             const double * __restrict __ATTR_ALIGN__(64) yre,
                                             const double * __restrict __ATTR_ALIGN__(64) yim,
                                             double * __restrict __ATTR_ALIGN__(64) zre,
                                             double * __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        rep  = _mm512_fmsub_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        zmm0 = _mm512_mul_pd(rep,rep);
                        zmm1 = _mm512_mul_pd(imp,imp);
                        zmm2 = _mm512_sqrt_pd(_mm512_add_pd(zmm0,zmm1));
                        _mm512_store_pd(&zre[0], _mm512_div_pd(rep,zmm2));
                        _mm512_store_pd(&zim[0], _mm512_div_pd(imp,zmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_zmm8r8(  const __m512d  xre,
                                             const __m512d  xim,
                                             const __m512d  yre,
                                             const __m512d  yim,
                                             __m512d * __restrict zre,
                                             __m512dd* __restrict zim) {

                        register __m512d rep,imp,zmm0,zmm1,zmm2;
                        rep  = _mm512_fmsub_pd(xre,yre,
                                               _mm512_mul_pd(xim,yim));
                        imp  = _mm512_fmadd_pd(xim,yre,
                                               _mm512_mul_pd(xre,yim));
                        zmm0 = _mm512_mul_pd(rep,rep);
                        zmm1 = _mm512_mul_pd(imp,imp);
                        zmm2 = _mm512_sqrt_pd(_mm512_add_pd(zmm0,zmm1));
                        *zre = _mm512_div_pd(rep,zmm2);
                        *zim = _mm512_div_pd(imp,zmm2);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_zmm8r8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm512_fmsub_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        sre  = _mm512_reduce_pd(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        sim  = _mm512_reduce_pd(imp);
                        *mim = sim*inv16;
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                             const double * __restrict __ATTR_ALIGN__(64) xim,
                                             const double * __restrict __ATTR_ALIGN__(64) yre,
                                             const double * __restrict __ATTR_ALIGN__(64) yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm512_fmsub_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        sre  = _mm512_reduce_pd(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        sim  = _mm512_reduce_pd(imp);
                        *mim = sim*inv16;
              } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_zmm8r8(const __m512d xre,
                                           const __m512d xim,
                                           const __m512d yre,
                                           const __m512d yim,
                                           double * __restrict mre,
                                           double * __restrict min) {

                        register __m512d rep,imp;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        sre = 0.0f;
                        rep  = _mm512_fmsub_pd(xre,yre,
                                               _mm512_mul_pd(xim,yim));
                        sre  = _mm512_reduce_pd(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm512_fmadd_pd(xim,yre,
                                               _mm512_mul_pd(xre,yim));
                        sim  = _mm512_reduce_pd(imp);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_zmm8r8_u(const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        sre  = 0.0;
                        rep  = _mm512_fmsub_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        sim  = 0.f;
                        den  = _mm512_fmadd_pd(zmm1,zmm1,
                                               _mm512_mul_pd(zmm3,zmm3));
                        rquot = _mm512_div_pd(rep,den);
                        sre   = _mm512_reduce_pd(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm512_div_pd(imp,den);
                        sim   = _mm512_reduce_pd(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                             const double * __restrict __ATTR_ALIGN__(64) xim,
                                             const double * __restrict __ATTR_ALIGN__(64) yre,
                                             const double * __restrict __ATTR_ALIGN__(64) yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        sre  = 0.0;
                        rep  = _mm512_fmsub_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        sim  = 0.0;
                        den  = _mm512_fmadd_pd(zmm1,zmm1,
                                               _mm512_mul_pd(zmm3,zmm3));
                        rquot = _mm512_div_pd(rep,den);
                        sre   = _mm512_reduce_pd(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm512_div_pd(imp,den);
                        sim   = _mm512_reduce_pd(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_zmm8r8(  const __m512d xre,
                                             const __m512d xim,
                                             const __m512d yre,
                                             const __m512d yim,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m512d rep,imp,den,rquot,iquot;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        sre  = 0.0;
                        rep  = _mm512_fmsub_pd(xre,yre,
                                               _mm512_mul_pd(xim,yim));
                        imp  = _mm512_fmadd_pd(xim,yre,
                                               _mm512_mul_pd(xre,yim));
                        sim  = 0.0;
                        den  = _mm512_fmadd_pd(yre,yre,
                                               _mm512_mul_pd(yim,yim));
                        rquot = _mm512_div_pd(rep,den);
                        sre   = _mm512_reduce_pd(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm512_div_pd(imp,den);
                        sim   = _mm512_reduce_pd(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_zmm8r8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim ) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp,magc1,magc2,vcmag;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        rep  = _mm512_fmadd_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        magc1= _mm512_mul_pd(rep,rep);
                        imp  = _mm512_fmsub_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        magc2= _mm512_mul_pd(imp,imp);
                        vcmag= _mm512_sqrt_pd(_mm512_add_pd(magc1,magc2));
                        _mm512_storeu_pd(&mre[0], _mm512_div_pd(rep,vcmag));
                        _mm512_storeu_pd(&mim[0], _mm512_div_pd(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                              const double * __restrict __ATTR_ALIGN__(64) xim,
                                              const double * __restrict __ATTR_ALIGN__(64) yre,
                                              const double * __restrict __ATTR_ALIGN__(64) yim,
                                              double * __restrict __ATTR_ALIGN__(64) mre,
                                              double * __restrict __ATTR_ALIGN__(64) mim ) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d rep,imp,magc1,magc2,vcmag;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        rep  = _mm512_fmad_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        magc1= _mm512_mul_pd(rep,rep);
                        imp  = _mm512_fmsub_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        magc2= _mm512_mul_pd(imp,imp);
                        vcmag= _mm512_sqrt_pd(_mm512_add_pd(magc1,magc2));
                        _mm512_store_pd(&mre[0], _mm512_div_pd(rep,vcmag));
                        _mm512_store_pd(&mim[0], _mm512_div_pd(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_zmm8r8(const __m512d xre,
                                            const __m512d xim,
                                            const __m512d yre,
                                            const __m512d yim,
                                            __m512d * __restrict mre,
                                            __m512d * __restrict mim) {

                        register __m512d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm512_fmad_pd(xre,yre,
                                               _mm512_mul_pd(xim,yim));
                        magc1= _mm512_mul_pd(rep,rep);
                        imp  = _mm512_fmsub_pd(xim,yre,
                                               _mm512_mul_pd(xre,yim));
                        magc2= _mm512_mul_pd(imp,imp);
                        vcmag= _mm512_sqrt_pd(_mm512_add_pd(magc1,magc2));
                        *mre = _mm512_div_pd(rep,vcmag);
                        *mim = _mm512_div_pd(imp,vcmag)
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_zmm8r8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim) {
                      
                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d re,im;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        re   = _mm512_fmadd_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        sre  = _mm512_reduce_pd(re);
                        *mre = sre*inv16;
                        im   = _mm512_fmsub_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        sim  = _mm512_reduce_pd(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) xre,
                                              const double * __restrict __ATTR_ALIGN__(64) xim,
                                              const double * __restrict __ATTR_ALIGN__(64) yre,
                                              const double * __restrict __ATTR_ALIGN__(64) yim,
                                              double * __restrict mre,
                                              double * __restrict mim) {
                      
                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d re,im;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        re   = _mm512_fmadd_pd(zmm0,zmm1,
                                               _mm512_mul_pd(zmm2,zmm3));
                        sre  = _mm512_reduce_pd(re);
                        *mre = sre*inv16;
                        im   = _mm512_fmsub_pd(zmm2,zmm1,
                                               _mm512_mul_pd(zmm0,zmm3));
                        sim  = _mm512_reduce_pd(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_zmm8r8(const __m512d xre,
                                            const __m512d xim,
                                            const __m512d yre,
                                            const __m512d yim,
                                            double * __restrict mre,
                                            double * __restrict min) {

                        register __m512d re,im;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        re   = _mm512_fmadd_pd(xre,yre,
                                               _mm512_mul_pd(xim,yim));
                        sre  = _mm512_reduce_pd(re);
                        *mre = sre*inv16;
                        im   = _mm512_fmsub_pd(xim,yre,
                                               _mm512_mul_pd(xre,yim));
                        sim  = _mm512_reduce_pd(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_zmm8r8_u(const double * __restrict xre,
                                              const double * __restrict xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        register __m512d re,im;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        re   = _mm512_loadu_pd(&xre[0]);
                        sre  = _mm512_reduce_pd(re);
                        *mre = sre*inv16;
                        im   = _mm512_loadu_pd(&xim[0]);
                        sim  = _mm512_reduce_pd(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_zmm8r8_u(const double * __restrict __ATTR_ALIGN__(64) xre,
                                              const double * __restrict __ATTR_ALIGN__(64) xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        register __m512d re,im;
                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        re   = _mm512_load_pd(&xre[0]);
                        sre  = _mm512_reduce_pd(re);
                        *mre = sre*inv16;
                        im   = _mm512_load_pd(&xim[0]);
                        sim  = _mm512_reduce_pd(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_zmm8r8(  const __m512d xre,
                                              const __m512d xim,
                                              double * __restrict mre,
                                              double * __restrict min) {

                        constexpr double inv16 = 0.0625;
                        double sre,sim;
                        sre  = _mm512_reduce_pd(xre);
                        *mre = sre*inv16;
                        sim  = _mm512_reduce_pd(xim);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_zmm8r8_u( const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre,
                                              double * __restrict mim ) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d re,im,cvmag;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        cvmag= _mm512_sqrt_pd(_mm512_fmadd_pd(zmm0,zmm1,
                                                              _mm512_mul_pd(zmm2,zmm3)));
                        _mm512_storeu_pd(&mre[0], _mm512_div_pd(zmm0,cvmag));
                        _mm512_storeu_pd(&mim[0], _mm512_div_pd(zmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) xre,
                                              const double * __restrict __ATTR_ALIGN__(64) xim,
                                              const double * __restrict __ATTR_ALIGN__(64) yre,
                                              const double * __restrict __ATTR_ALIGN__(64) yim,
                                              double * __restrict __ATTR_ALIGN__(64) mre,
                                              double * __restrict __ATTR_ALIGN__(64) mim ) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d re,im,cvmag;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        cvmag= _mm512_sqrt_pd(_mm512_fmadd_pd(zmm0,zmm1,
                                                              _mm512_mul_pd(zmm2,zmm3)));
                        _mm512_store_pd(&mre[0], _mm512_div_pd(zmm0,cvmag));
                        _mm512_store_pd(&mim[0], _mm512_div_pd(zmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_zmm8r8( const __m512d xre,
                                            const __m512d xim,
                                            const __m512d yre,
                                            const __m512d yim,
                                            __m512d * __restrict mre,
                                            __m512d * __restrict mim ) {

                        register __m512d re,im,cvmag;
                        cvmag= _mm512_sqrt_pd(_mm512_fmadd_pd(xre,yre,
                                                    _mm512_mul_pd(xim,yim)));
                        *mre = _mm512_div_pd(xre,cvmag));
                        *mim =  _mm512_div_pd(xim,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_zmm8r8_u( const double * __restrict xre,
                                              const double * __restrict xim,
                                              const double * __restrict yre,
                                              const double * __restrict yim,
                                              double * __restrict mre) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d cvmag;
                        zmm0 = _mm512_loadu_pd(&xre[0]);
                        zmm1 = _mm512_loadu_pd(&yre[0]);
                        zmm2 = _mm512_loadu_pd(&xim[0]);
                        zmm3 = _mm512_loadu_pd(&yim[0]);
                        cvmag= _mm512_sqrt_pd(_mm512_fmadd_pd(zmm0,zmm1,
                                                          _mm512_mul_pd(zmm2,zmm3)));
                        _mm512_storeu_pd(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) xre,
                                              const double * __restrict __ATTR_ALIGN__(64) xim,
                                              const double * __restrict __ATTR_ALIGN__(64) yre,
                                              const double * __restrict __ATTR_ALIGN__(64) yim,
                                              double * __restrict __ATTR_ALIGN__(64) mre) {

                        register __m512d zmm0,zmm1,zmm2,zmm3;
                        register __m512d cvmag;
                        zmm0 = _mm512_load_pd(&xre[0]);
                        zmm1 = _mm512_load_pd(&yre[0]);
                        zmm2 = _mm512_load_pd(&xim[0]);
                        zmm3 = _mm512_load_pd(&yim[0]);
                        cvmag= _mm512_sqrt_pd(_mm512_fmadd_pd(zmm0,zmm1,
                                                          _mm512_mul_pd(zmm2,zmm3)));
                        _mm512_store_pd(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_zmm8r8(   const __m512d xre,
                                              const __m512d xim,
                                              const __m512d yre,
                                              const __m512d yim,
                                              __m512d * __restrict  mre) {

                        register __m512d cvmag;
                        cvmag= _mm512_sqrt_pd(_mm512_fmadd_pd(xre,yre,
                                                          _mm512_mul_pd(xim,yim)));
                        *mre = cvmag;
             }






                                         


      } // math


} // gms















#endif /*__GMS_COMPLEX_ZMM8R8_HPP__*/
