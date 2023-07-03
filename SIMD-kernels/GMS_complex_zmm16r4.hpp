
#ifndef __GMS_COMPLEX_ZMM16R4_HPP__
#define __GMS_COMPLEX_ZMM16R4_HPP__ 251220220954


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

    const unsigned int GMS_COMPLEX_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_COMPLEX_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_COMPLEX_ZMM16R4_FULLVER =
      1000U*GMS_COMPLEX_ZMM16R4_MAJOR+
      100U*GMS_COMPLEX_ZMM16R4_MINOR+
      10U*GMS_COMPLEX_ZMM16R4_MICRO;
    const char * const GMS_COMPLEX_ZMM16R4_CREATION_DATE = "25-12-2022 09:54 AM +00200 (SUN 25 DEC 2022 GMT+2)";
    const char * const GMS_COMPLEX_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex number implementation."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace  gms {


       namespace math {
       
       
                   struct __ATTR_ALIGN__(64) zmm16c4_t {
                   
                          __m512 re;
                          __m512 im;
                   };

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_ps(&xre[0]);
                        zmm1  = _mm512_loadu_ps(&yre[0]);
                        _mm512_storeu_ps(&zre[0], _mm512_add_ps(zmm0,zmm1));
                        zmm2  = _mm512_loadu_ps(&xim[0]);
                        zmm3  = _mm512_loadu_ps(&yim[0]);
                        _mm512_storeu_ps(&zim[0], _mm512_add_ps(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       const float * __restrict __ATTR_ALIGN__(64) yre,
                                       const float * __restrict __ATTR_ALIGN__(64) yim,
                                       float *       __restrict __ATTR_ALIGN__(64) zre,
                                       float *       __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_ps(&xre[0]);
                        zmm1  = _mm512_load_ps(&yre[0]);
                        _mm512_store_ps(&zre[0], _mm512_add_ps(zmm0,zmm1));
                        zmm2  = _mm512_load_ps(&xim[0]);
                        zmm3  = _mm512_load_ps(&yim[0]);
                        _mm512_store_ps(&zim[0], _mm512_add_ps(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 yre,
                                     const __m512 yim,
                                     __m512 * __restrict zre,
                                     __m512 * __restrict zim) {
                     
                        register __m512 zmm0,zmm1;
                        zmm0  = _mm512_add_ps(xre,yre);
                        *zre  = zmm0;
                        zmm1  = _mm512_add_ps(xim,yim);
                        *zim  = zmm1;
                }
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cadd_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 yre,
                                          const __m512 yim) {
                                     
                        zmm16c4_t cv;
                        register __m512 zmm0,zmm1;
                        zmm0   = _mm512_add_ps(xre,yre);
                        cv.re  = zmm0;
                        zmm1   = _mm512_add_ps(xim,yim);
                        cv.im  = zmm1;  
                        return (cv);            
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 s,
                                     __m512 * __restrict     zre,
                                     __m512 * __restrict     zim) {

                        *zre = _mm512_add_ps(xre,s);
                        *zim = _mm512_add_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           zmm16c4_t cadd_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 s) {
                      
                      zmm16c4_t cv;
                      cv.re =  _mm512_add_ps(xre,s);
                      cv.im =  _mm512_add_ps(xim,s);
                      return (cv);                       
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm16r4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_ps(&xre[0]);
                        zmm1  = _mm512_loadu_ps(&xim[0]);
                        zmm2  = _mm512_loadu_ps(&zre[0]);
                        zmm3  = _mm512_loadu_ps(&zim[0])
                        _mm512_storeu_ps(&zre[0], _mm512_add_ps(zmm2,zmm0));
                        _mm512_storeu_ps(&zim[0], _mm512_add_ps(zmm3,zmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cadd_zmm16r4_aip(const float * __restrict __ATTR_ALIGN__(64) xre,
                                         const float * __restrict __ATTR_ALIGN__(64) xim,
                                         float *       __restrict __ATTR_ALIGN__(64) zre,
                                         float *       __restrict __ATTR_ALIGN__(64) zim) {
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_ps(&xre[0]);
                        zmm1  = _mm512_load_ps(&xim[0]);
                        zmm2  = _mm512_load_ps(&zre[0]);
                        zmm3  = _mm512_load_ps(&zim[0])
                        _mm512_store_ps(&zre[0], _mm512_add_ps(zmm2,zmm0));
                        _mm512_store_ps(&zim[0], _mm512_add_ps(zmm3,zmm1));
              }


                ////////////////////////////////////////////////////////////////////

                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_ps(&xre[0]);
                        zmm1  = _mm512_loadu_ps(&yre[0]);
                        _mm512_storeu_ps(&zre[0], _mm512_sub_ps(zmm0,zmm1));
                        zmm2  = _mm512_loadu_ps(&xim[0]);
                        zmm3  = _mm512_loadu_ps(&yim[0]);
                        _mm512_storeu_ps(&zim[0], _mm512_sub_ps(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       const float * __restrict __ATTR_ALIGN__(64) yre,
                                       const float * __restrict __ATTR_ALIGN__(64) yim,
                                       float *       __restrict __ATTR_ALIGN__(64) zre,
                                       float *       __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_ps(&xre[0]);
                        zmm1  = _mm512_load_ps(&yre[0]);
                        _mm512_store_ps(&zre[0], _mm512_sub_ps(zmm0,zmm1));
                        zmm2  = _mm512_load_ps(&xim[0]);
                        zmm3  = _mm512_load_ps(&yim[0]);
                        _mm512_store_ps(&zim[0], _mm512_sub_ps(zmm2,zmm3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 yre,
                                     const __m512 yim,
                                     __m512 * __restrict     zre,
                                     __m512 * __restrict     zim) {
                     
                        register __m512 zmm0,zmm1;
                        zmm0  = _mm512_sub_ps(xre,yre);
                        *zre  = zmm0;
                        zmm1  = _mm512_sub_ps(xim,yim);
                        *zim  = zmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t csub_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 yre,
                                          const __m512 yim) {
                                    
                        zmm16c4_t cv;
                        register __m512 zmm0,zmm1;
                        zmm0  = _mm512_sub_ps(xre,yre);
                        cv.re  = zmm0;
                        zmm1  = _mm512_sub_ps(xim,yim);
                        cv.im  = zmm1;
                        return (cv);
                }
                
                


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 s,
                                     __m512 * __restrict     zre,
                                     __m512 * __restrict     zim) {

                        *zre = _mm512_sub_ps(xre,s);
                        *zim = _mm512_sub_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t csub_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 s) {
                                    
                        zmm16c4_t cv;
                        cv.re = _mm512_sub_ps(xre,s);
                        cv.im = _mm512_sub_ps(xim,s);
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm16r4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_ps(&xre[0]);
                        zmm1  = _mm512_loadu_ps(&xim[0]);
                        zmm2  = _mm512_loadu_ps(&zre[0]);
                        zmm3  = _mm512_loadu_ps(&zim[0])
                        _mm512_storeu_ps(&zre[0], _mm512_sub_ps(zmm2,zmm0));
                        _mm512_storeu_ps(&zim[0], _mm512_sub_ps(zmm3,zmm1));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csub_zmm16r4_aip(const float * __restrict __ATTR_ALIGN__(64) xre,
                                         const float * __restrict __ATTR_ALIGN__(64) xim,
                                         float *       __restrict __ATTR_ALIGN__(64) zre,
                                         float *       __restrict __ATTR_ALIGN__(64) zim) {
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_ps(&xre[0]);
                        zmm1  = _mm512_load_ps(&xim[0]);
                        zmm2  = _mm512_load_ps(&zre[0]);
                        zmm3  = _mm512_load_ps(&zim[0])
                        _mm512_store_ps(&zre[0], _mm512_sub_ps(zmm2,zmm0));
                        _mm512_store_ps(&zim[0], _mm512_sub_ps(zmm3,zmm1));
              }


               ////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                           register __m512 zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_loadu_ps(&xre[0]);
                           zmm1  = _mm512_loadu_ps(&yre[0]);
                           zmm2  = _mm512_loadu_ps(&xim[0]);
                           zmm3  = _mm512_loadu_ps(&yim[0]);
                           zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                           _mm512_storeu_ps(&zre[0], zmm4);
                           zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                           _mm512_storeu_ps(&zim[0], zmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       const float * __restrict __ATTR_ALIGN__(64) yre,
                                       const float * __restrict __ATTR_ALIGN__(64) yim,
                                       float *       __restrict __ATTR_ALIGN__(64) zre,
                                       float *       __restrict __ATTR_ALIGN__(64) zim) {

                           register __m512 zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_load_ps(&xre[0]);
                           zmm1  = _mm512_load_ps(&yre[0]);
                           zmm2  = _mm512_load_ps(&xim[0]);
                           zmm3  = _mm512_load_ps(&yim[0]);
                           zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                                        _mm512_mul_ps(zmm2,zmm3));
                           _mm512_store_ps(&zre[0], zmm4);
                           zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                                        _mm512_mul_ps(zmm0,zmm3));
                           _mm512_store_ps(&zim[0], zmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 yre,
                                     const __m512 yim,
                                     __m512 * __restrict     zre,
                                     __m512 * __restrict     zim) {

                         register __m512 zmm0,zmm1;
                         zmm0 = _mm512_sub_ps(_mm512_mul_ps(xre,yre),
                                              _mm512_mul_ps(xim,yim));
                         *zre  = zmm0;
                         zmm1 = _mm512_mul_ps(_mm512_mul_ps(xim,yre),
                                              _mm512_mul_ps(xre,yim));
                         *zim  = zmm1;
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cmul_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 yre,
                                          const __m512 yim) {
                                     
                         zmm16c4_t cv
                         register __m512 zmm0,zmm1;
                         zmm0 = _mm512_sub_ps(_mm512_mul_ps(xre,yre),
                                              _mm512_mul_ps(xim,yim));
                         cv.re  = zmm0;
                         zmm1 = _mm512_mul_ps(_mm512_mul_ps(xim,yre),
                                              _mm512_mul_ps(xre,yim));
                         cv.im  = zmm1;
                         return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 s,
                                     __m512 * __restrict   zre,
                                     __m512 * __restrict   zim) {

                        *zre = _mm512_mul_ps(xre,s);
                        *zim = _mm512_mul_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cmul_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 s) {
                                     
                        zmm16c4_t cv;
                        cv.re = _mm512_mul_ps(xre,s);
                        cv.im = _mm512_mul_ps(xim,s);
                        return (cv);
               }
               


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm16r4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {

                           register __m512 zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_loadu_ps(&xre[0]);
                           zmm1  = _mm512_loadu_ps(&zre[0]);
                           zmm2  = _mm512_loadu_ps(&xim[0]);
                           zmm3  = _mm512_loadu_ps(&zim[0]);
                           zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                 _mm512_mul_ps(zmm2,zmm3));
                           _mm512_storeu_ps(&zre[0], zmm4);
                           zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                 _mm512_mul_ps(zmm0,zmm3));
                           _mm512_storeu_ps(&zim[0], zmm5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmul_zmm16r4_aip(const float * __restrict __ATTR_ALIGN__(64) xre,
                                         const float * __restrict __ATTR_ALIGN__(64) xim,
                                         float *       __restrict __ATTR_ALIGN__(64) zre,
                                         float *       __restrict __ATTR_ALIGN__(64) zim) {

                           register __m512 zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
                           zmm0  = _mm512_load_ps(&xre[0]);
                           zmm1  = _mm512_load_ps(&zre[0]);
                           zmm2  = _mm512_load_ps(&xim[0]);
                           zmm3  = _mm512_load_ps(&zim[0]);
                           zmm4  = _mm512_sub_ps(_mm512_mul_ps(zmm0,zmm1),
                                                 _mm512_mul_ps(zmm2,zmm3));
                           _mm512_store_ps(&zre[0], zmm4);
                           zmm5  = _mm512_mul_ps(_mm512_mul_ps(zmm2,zmm1),
                                                 _mm512_mul_ps(zmm0,zmm3));
                           _mm512_store_ps(&zim[0], zmm5);
               }

                 ////////////////////////////////////////////////////////////////////////////

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3; 
                        register __m512 zmm4,zmm5,zmm6;
                        zmm0  = _mm512_loadu_ps(&xre[0]); //a
                        zmm1  = _mm512_loadu_ps(&yim[0]); //d
                        zmm2  = _mm512_loadu_ps(&xim[0]); //b
                        zmm3  = _mm512_loadu_ps(&yre[0]); //c
                        zmm4  = _mm512_fmadd_ps(zmm0,zmm3,
                                                _mm512_mul_ps(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_ps(zmm2,zmm3,
                                                _mm512_mul_ps(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_ps(zmm3,zmm3),
                                                _mm512_mul_ps(zmm1,zmm1));
                        _mm512_storeu_ps(&zre[0], _mm512_div_ps(zmm4,zmm6));
                        _mm512_storeu_ps(&zim[0], _mm512_div_ps(zmm5,zmm6));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4_a(const float * __restrict xre,
                                       const float * __restrict xim,
                                       const float * __restrict yre,
                                       const float * __restrict yim,
                                       float *       __restrict zre,
                                       float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3; 
                        register __m512 zmm4,zmm5,zmm6;
                        zmm0  = _mm512_load_ps(&xre[0]); //a
                        zmm1  = _mm512_load_ps(&yim[0]); //d
                        zmm2  = _mm512_load_ps(&xim[0]); //b
                        zmm3  = _mm512_load_ps(&yre[0]); //c
                        zmm4  = _mm512_fmadd_ps(zmm0,zmm3,
                                                _mm512_mul_ps(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_ps(zmm2,zmm3,
                                                _mm512_mul_ps(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_ps(zmm3,zmm3,
                                                _mm512_mul_ps(zmm1,zmm1));
                        _mm512_store_ps(&zre[0], _mm512_div_ps(zmm4,zmm6));
                        _mm512_store_ps(&zim[0], _mm512_div_ps(zmm5,zmm6));
                }
              
                                         
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 yre,
                                     const __m512 yim,
                                     __m512 * __restrict zre,
                                     __m512 * __restrict zim) {

                      register __m512 zmm0,zmm1,zmm2;
                      zmm0 = _mm512_fmadd_ps(xre,yre,
                                           _mm512_mul_ps(xim,yim));
                      zmm1 = _mm512_fmsub_ps(xim,yre,
                                           _mm512_mul_ps(xre,yim));
                      zmm2 = _mm512_fmadd_ps(zmm3,zmm3,
                                           _mm512_mul_ps(zmm1,zmm1));
                      *zre  = _mm512_div_ps(zmm0,zmm2);
                      *zim  = _mm512_div_ps(zmm1,zmm2);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cdiv_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 yre,
                                          const __m512 yim) {
                                     
                      zmm16c4_t
                      register __m512 zmm0,zmm1,zmm2;
                      zmm0 = _mm512_fmadd_ps(xre,yre,
                                           _mm512_mul_ps(xim,yim));
                      zmm1 = _mm512_fmsub_ps(xim,yre,
                                           _mm512_mul_ps(xre,yim));
                      zmm2 = _mm512_fmadd_ps(zmm3,zmm3,
                                           _mm512_mul_ps(zmm1,zmm1));
                      cv.re  = _mm512_div_ps(zmm0,zmm2);
                      cv.im  = _mm512_div_ps(zmm1,zmm2);
                      return (cv);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     const __m512 s,
                                     __m512 * __restrict zre,
                                     __m512 * __restrict zim) {

                        *zre = _mm512_div_ps(xre,s);
                        *zim = _mm512_div_ps(xim,s);
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cdiv_zmm16r4(const __m512 xre,
                                          const __m512 xim,
                                          const __m512 s) {
                                     
                         zmm16c4_t cv;
                         cv.re = _mm512_div_ps(xre,s);
                         cv.im = _mm512_div_ps(xim,s);
                         return (cv);
               }
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4_s(const __m512 s,
                                       const __m512 xre,
                                       const __m512 xim,
                                       __m512 * __restrict zre,
                                       __m512 * __restrict zim) {
                        
                        register __m512 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_ps();
                        cdiv_zmm16r4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cdiv_zmm16r4_s(const __m512 s,
                                            const __m512 xre,
                                            const __m512 xim) {
                                       
                        zmm16c4_t cv;
                        register __m512 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_ps();
                        cdiv_zmm16r4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re  = tmpr;
                        cv.zim = tmpi;
                        return (cv);                     
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm16r4_s(const __m512 s,
                                             const __m512 xre,
                                             const __m512 xim,
                                             __m512 * __restrict zre,
                                             __m512 * __restrict zim) {
                                             
                        register __m512 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_ps(); 
                        cdiv_smith_zmm16r4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        *zre = tmpr;
                        *zim = tmpi;                   
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cdiv_smith_zmm16r4_s(const __m512 s,
                                                  const __m512 xre,
                                                  const __m512 xim) {
                                             
                        zmm16c4_t cv;                    
                        register __m512 t0r,t0i,tmpr,tmpi;
                        t0r = s;
                        t0i = _mm512_setzero_ps(); 
                        cdiv_smith_zmm16r4(t0r,t0i,xre,xim,&tmpr,&tmpi);
                        cv.re = tmpr;
                        cv.im = tmpi;  
                        return (cv);                 
                 }
                 
                 
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4_uip(const float * __restrict xre,
                                         const float * __restrict xim,
                                         float *       __restrict zre,
                                         float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3; 
                        register __m512 zmm4,zmm5,zmm6;
                        zmm0  = _mm512_loadu_ps(&xre[0]); //a
                        zmm1  = _mm512_loadu_ps(&zim[0]); //d
                        zmm2  = _mm512_loadu_ps(&xim[0]); //b
                        zmm3  = _mm512_loadu_ps(&zre[0]); //c
                        zmm4  = _mm512_fmadd_ps(zmm0,zmm3,
                                                _mm512_mul_ps(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_ps(zmm2,zmm3,
                                                _mm512_mul_ps(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_ps(zmm3,zmm3,
                                                _mm512_mul_ps(zmm1,zmm1));
                        _mm512_storeu_ps(&zre[0], _mm512_div_ps(zmm4,zmm6));
                        _mm512_storeu_ps(&zim[0], _mm512_div_ps(zmm5,zmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_zmm16r4_aip(const float * __restrict __ATTR_ALIGN__(64) xre,
                                         const float * __restrict __ATTR_ALIGN__(64) xim,
                                         float *       __restrict __ATTR_ALIGN__(64) zre,
                                         float *       __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3; 
                        register __m512 zmm4,zmm5,zmm6;
                        zmm0  = _mm512_load_ps(&xre[0]); //a
                        zmm1  = _mm512_load_ps(&zim[0]); //d
                        zmm2  = _mm512_load_ps(&xim[0]); //b
                        zmm3  = _mm512_load_ps(&zre[0]); //c
                        zmm4  = _mm512_fmadd_ps(zmm0,zmm3,
                                                _mm512_mul_ps(zmm2,zmm1));
                        zmm5  = _mm512_fmsub_ps(zmm2,zmm3,
                                                _mm512_mul_ps(zmm0,zmm1));
                        zmm6  = _mm512_fmadd_ps(zmm3,zmm3,
                                                _mm512_mul_ps(zmm1,zmm1));
                        _mm512_store_ps(&zre[0], _mm512_div_ps(zmm4,zmm6));
                        _mm512_store_ps(&zim[0], _mm512_div_ps(zmm5,zmm6));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm16r4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float *       __restrict zre,
                                             float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 r,den;
                        __mmask16 m = 0x0;
                        zmm0 = _mm512_loadu_ps(&yre[0]); // c
                        zmm1 = _mm512_loadu_ps(&yim[0]); // d
                        zmm2 = _mm512_loadu_ps(&xre[0]); // a
                        zmm3 = _mm512_loadu_ps(&xim[0]); // b
                        m    = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                  _mm512_abs_ps(zmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_ps(m,_mm512_div_ps(zmm0,zmm1),
                                                      _mm512_div_ps(zmm1,zmm0)); // r
                        den  = _mm512_mask_blend_ps(m,_mm512_fmadd_ps(r,zmm0,zmm1),
                                                      _mm512_fmadd_ps(r,zmm1,zmm0));
                        _mm512_storeu_ps(&zre[0], _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmadd_ps(zmm2,r,zmm3),den),
                                                _mm512_div_ps(_mm512_fmadd_ps(zmm3,r,zmm2),den)));
                        _mm512_storeu_ps(&zim[0], _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmsub_ps(zmm3,r,zmm2),den),
                                                _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm16r4_a(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float *       __restrict zre,
                                             float *       __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 r,den;
                        __mmask16 m = 0x0;
                        zmm0 = _mm512_load_ps(&yre[0]); // c
                        zmm1 = _mm512_load_ps(&yim[0]); // d
                        zmm2 = _mm512_load_ps(&xre[0]); // a
                        zmm3 = _mm512_load_ps(&xim[0]); // b
                        m    = _mm512_cmp_ps_mask(_mm512_abs_ps(zmm0),
                                                  _mm512_abs_ps(zmm1),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_ps(m,_mm512_div_ps(zmm0,zmm1),
                                                      _mm512_div_ps(zmm1,zmm0)); // r
                        den  = _mm512_mask_blend_ps(m,_mm512_fmadd_ps(r,zmm0,zmm1),
                                                      _mm512_fmadd_ps(r,zmm1,zmm0));
                        _mm512_storeu_ps(&zre[0], _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmadd_ps(zmm2,r,zmm3),den),
                                                _mm512_div_ps(_mm512_fmadd_ps(zmm3,r,zmm2),den)));
                        _mm512_storeu_ps(&zim[0], _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmsub_ps(zmm3,r,zmm2),den),
                                                _mm512_div_ps(_mm512_sub_ps(zmm3,_mm512_mul_ps(zmm2,r)),den)));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cdiv_smith_zmm16r4(const __m512 xre,
                                           const __m512 xim,
                                           const __m512 yre,
                                           const __m512 yim,
                                           __m512 * __restrict zre,
                                           __m512 * __restrict zim) {

                        register __m512 r,den;
                        __mmask16 m = 0x0;
                        m    = _mm512_cmp_ps_mask(_mm512_abs_ps(yre),
                                                  _mm512_abs_ps(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_ps(m,_mm512_div_ps(yre,yim),
                                                      _mm512_div_ps(yim,yre)); // r
                        den  = _mm512_mask_blend_ps(m,_mm512_fmadd_ps(r,yre,yim),
                                                      _mm512_fmadd_ps(r,yim,yre));
                        *zre  =  _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmadd_ps(xre,r,xim),den),
                                                _mm512_div_ps(_mm512_fmadd_ps(xim,r,xre),den));
                        *zim  =  _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmsub_ps(xim,r,xre),den),
                                                _mm512_div_ps(_mm512_sub_ps(xim,_mm512_mul_ps(xre,r)),den)));
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cdiv_smith_zmm16r4(const __m512 xre,
                                                const __m512 xim,
                                                const __m512 yre,
                                                const __m512 yim) {
                                           
                        zmm16c4_t cv
                        register __m512 r,den;
                        __mmask16 m = 0x0;
                        m    = _mm512_cmp_ps_mask(_mm512_abs_ps(yre),
                                                  _mm512_abs_ps(yim),
                                                  _CMP_GE_OQ);
                        r    = _mm512_mask_blend_ps(m,_mm512_div_ps(yre,yim),
                                                      _mm512_div_ps(yim,yre)); // r
                        den  = _mm512_mask_blend_ps(m,_mm512_fmadd_ps(r,yre,yim),
                                                      _mm512_fmadd_ps(r,yim,yre));
                        cv.re  =  _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmadd_ps(xre,r,xim),den),
                                                _mm512_div_ps(_mm512_fmadd_ps(xim,r,xre),den));
                        cv.im  =  _mm512_mask_blend_ps(m,
                                                _mm512_div_ps(_mm512_fmsub_ps(xim,r,xre),den),
                                                _mm512_div_ps(_mm512_sub_ps(xim,_mm512_mul_ps(xre,r)),den)));
                        return (cv);
               }


#include "GMS_sleefsimdsp.hpp"


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_zmm16r4_u(const float * __restrict re,
                                       const float * __restrict im,
                                       float * __restrict  cabs) {

                        register __m512 zmm0,zmm1,zmm2,zmm3,zmm4;
                        zmm0  = _mm512_loadu_ps(&re[0]);
                        zmm1  = _mm512_mul_ps(zmm0,zmm0);
                        zmm2  = _mm512_loadu_ps(&im[0]);
                        zmm3  = _mm512_mul_ps(zmm2,zmm2);
                        zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                        _mm512_storeu_ps(&cabs[0],zmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cabs_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                       const float * __restrict __ATTR_ALIGN__(64) im,
                                       float * __restrict  __ATTR_ALIGN__(64) cabs) {

                        register __m512 zmm0,zmm1,zmm2,zmm3,zmm4;
                        zmm0  = _mm512_load_ps(&re[0]);
                        zmm1  = _mm512_mul_ps(zmm0,zmm0);
                        zmm2  = _mm512_load_ps(&im[0]);
                        zmm3  = _mm512_mul_ps(zmm2,zmm2);
                        zmm4  = xsqrtf(_mm512_add_ps(zmm1,zmm3));
                        _mm512_store_ps(&cabs[0],zmm4);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 cabs_zmm16r4(const __m512 re,
                                       const __m512 im) {

                        register __m512 zmm0,zmm1,cabs;
                        zmm0 = _mm512_mul_ps(re,re);
                        zmm1 = _mm512_mul_ps(im,im);
                        cabs = xsqrtf(_mm512_add_ps(zmm0,zmm1));
                        return (cabs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_zmm16r4_u(const float * __restrict re,
                                       const float * __restrict im,
                                       float * __restrict  carg) {

                        register __m512 zmm0,zmm1;
                        zmm0 = _mm512_loadu_ps(&re[0]);
                        zmm1 = _mm512_loadu_ps(&im[0]);
                        _mm512_storeu_ps(&carg[0], xatan2f(zmm0,zmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void carg_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                       const float * __restrict __ATTR_ALIGN__(64) im,
                                       float * __restrict  __ATTR_ALIGN__(64) carg) {

                        register __m512 zmm0,zmm1;
                        zmm0 = _mm512_load_ps(&re[0]);
                        zmm1 = _mm512_load_ps(&im[0]);
                        _mm512_store_ps(&carg[0], xatan2f(zmm0,zmm1));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 carg_zmm16r4(const __m512 re,
                                       const __m512 im) {

                       register __m512 carg;
                       carg = xatan2f(re,im);
                       return (carg);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm16r4_u(float * __restrict re,
                                        float * __restrict im) {

                        register __m512 c;
                        c = negate_zmm16r4(_mm512_loadu_ps(&im[0]));
                        _mm512_storeu_ps(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) re,
                                        float * __restrict __ATTR_ALIGN__(64) im) {
                                        
                        register __m512 c;
                        c = negate_zmm16r4(_mm512_load_ps(&im[0]));
                        _mm512_store_ps(&im[0],c);
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm16r4(__m512 * __restrict re,
                                      __m512 * __restrict im) {
                         
                        register __m512 c;              
                        c = negate_zmm16r4(*im);
                        *im = c;
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm16r4_v2(const __m512 xre,
                                         const __m512 xim,
                                         __m512 * __restrict yre,
                                         __m512 * __restrict yim) {
                         
                        //register __m512 c;              
                        //c = negate_zmm16r4(*im);
                        //*im = c;
                        *yre = xre; 
                        *yim = negate_zmm16r4(xim);
                   } 
                   
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cconj_zmm16r4_v2(const __m512 xre,
                                              const __m512 xim) {                                              
                         
                        //register __m512 c;              
                        //c = negate_zmm16r4(*im);
                        //*im = c;
                        zmm16c4_t cv;
                        cv.re = xre; 
                        cv.im = negate_zmm16r4(xim);
                        return (cv);
                   } 
                   
                   


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict  csre,
                                       float * __restrict  csim) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_loadu_ps(&xre[0]);
                      zmm1  = _mm512_loadu_ps(&xim[0]);
                      zmm2  = _mm512_mul_ps(xcosf(zmm0),xcoshf(zmm1));
                      _mm512_storeu_ps(&csre[0],zmm2);
                      zmm3  = _mm512_mul_ps(xsinf(zmm0),xsinhf(zmm1));
                      _mm512_storeu_ps(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       float * __restrict  __ATTR_ALIGN__(64) csre,
                                       float * __restrict  __ATTR_ALIGN__(64) csim) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_load_ps(&xre[0]);
                      zmm1  = _mm512_load_ps(&xim[0]);
                      zmm2  = _mm512_mul_ps(xcosf(zmm0),xcoshf(zmm1));
                      _mm512_store_ps(&csre[0],zmm2);
                      zmm3  = _mm512_mul_ps(xsinf(zmm0),xsinhf(zmm1));
                      _mm512_store_ps(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccos_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     __m512 * __restrict csre,
                                     __m512 * __restrict csim) {

                      register __m512 zmm0,zmm1;
                      zmm0  = _mm512_mul_ps(xcosf(xre),xcoshf(xim));
                      *csre = zmm0;
                      zmm1  = _mm512_mul_ps(xsinf(xre),xsinhf(xim));
                      *csim = zmm1; 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t ccos_zmm16r4(const __m512 xre,
                                          const __m512 xim) {
                                    
                      zmm16c4_t cv;
                      register __m512 zmm0,zmm1;
                      zmm0  = _mm512_mul_ps(xcosf(xre),xcoshf(xim));
                      cv.re = zmm0;
                      zmm1  = _mm512_mul_ps(xsinf(xre),xsinhf(xim));
                      cv.im = zmm1;
                      return (cv); 
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict  csre,
                                       float * __restrict  csim) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_loadu_ps(&xre[0]);
                      zmm1  = _mm512_loadu_ps(&xim[0]);
                      zmm2  = _mm512_mul_ps(xcoshf(zmm0),xcosf(zmm1));
                      _mm512_storeu_ps(&csre[0],zmm2);
                      zmm3  = _mm512_mul_ps(xsinhf(zmm0),xsinf(zmm1));
                      _mm512_storeu_ps(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       float * __restrict  __ATTR_ALIGN__(64) csre,
                                       float * __restrict  __ATTR_ALIGN__(64) csim) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0  = _mm512_load_ps(&xre[0]);
                      zmm1  = _mm512_load_ps(&xim[0]);
                      zmm2  = _mm512_mul_ps(xcoshf(zmm0),xcosf(zmm1));
                      _mm512_store_ps(&csre[0],zmm2);
                      zmm3  = _mm512_mul_ps(xsinhf(zmm0),xsinf(zmm1));
                      _mm512_store_ps(&csim[0],zmm3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ccosh_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     __m512 * __restrict csre,
                                     __m512 * __restrict csim) {

                      register __m512 zmm0,zmm1;
                      zmm0  = _mm512_mul_ps(xcoshf(xre),xcosf(xim));
                      *csre = zmm0;
                      zmm1  = _mm512_mul_ps(xsinhf(xre),xsinf(xim));
                      *csim = zmm1; 
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t ccosh_zmm16r4(const __m512 xre,
                                           const __m512 xim) {
                                          
                      zmm16c4_t cv;
                      register __m512 zmm0,zmm1;
                      zmm0  = _mm512_mul_ps(xcoshf(xre),xcosf(xim));
                      cv.re = zmm0;
                      zmm1  = _mm512_mul_ps(xsinhf(xre),xsinf(xim));
                      cv.im = zmm1; 
                      return (cv);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_zmm16r4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_ps(&xre[0]);
                      zmm1 = _mm512_loadu_ps(&yre[0]);
                      _mm512_storeu_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                      zmm2 = _mm512_loadu_ps(&xim[0]);
                      zmm3 = _mm512_loadu_ps(&yim[0]);
                      _mm512_storeu_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                      const float * __restrict __ATTR_ALIGN__(64) xim,
                                      const float * __restrict __ATTR_ALIGN__(64) yre,
                                      const float * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_ps(&xre[0]);
                      zmm1 = _mm512_load_ps(&yre[0]);
                      _mm512_store_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_EQ_OQ));
                      zmm2 = _mm512_load_ps(&xim[0]);
                      zmm3 = _mm512_load_ps(&yim[0]);
                      _mm512_store_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_EQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ceq_zmm16r4(const __m512 xre,
                                    const __m512 xim,
                                    const __m512 yre,
                                    const __m512 yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_ps_mask(xre,yre,_CMP_EQ_OQ);
                         *eqi = _mm512_cmp_ps_mask(xim,yim,_CMP_EQ_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_zmm16r4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_ps(&xre[0]);
                      zmm1 = _mm512_loadu_ps(&yre[0]);
                      _mm512_storeu_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_GT_OQ));
                      zmm2 = _mm512_loadu_ps(&xim[0]);
                      zmm3 = _mm512_loadu_ps(&yim[0]);
                      _mm512_storeu_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_GT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                      const float * __restrict __ATTR_ALIGN__(64) xim,
                                      const float * __restrict __ATTR_ALIGN__(64) yre,
                                      const float * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_ps(&xre[0]);
                      zmm1 = _mm512_load_ps(&yre[0]);
                      _mm512_store_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_GT_OQ));
                      zmm2 = _mm512_load_ps(&xim[0]);
                      zmm3 = _mm512_load_ps(&yim[0]);
                      _mm512_store_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_GT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cgt_zmm16r4(const __m512 xre,
                                    const __m512 xim,
                                    const __m512 yre,
                                    const __m512 yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_ps_mask(xre,yre,_CMP_GT_OQ);
                         *eqi = _mm512_cmp_ps_mask(xim,yim,_CMP_GT_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_zmm16r4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_ps(&xre[0]);
                      zmm1 = _mm512_loadu_ps(&yre[0]);
                      _mm512_storeu_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_LT_OQ));
                      zmm2 = _mm512_loadu_ps(&xim[0]);
                      zmm3 = _mm512_loadu_ps(&yim[0]);
                      _mm512_storeu_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_LT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                      const float * __restrict __ATTR_ALIGN__(64) xim,
                                      const float * __restrict __ATTR_ALIGN__(64) yre,
                                      const float * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_ps(&xre[0]);
                      zmm1 = _mm512_load_ps(&yre[0]);
                      _mm512_store_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_LT_OQ));
                      zmm2 = _mm512_load_ps(&xim[0]);
                      zmm3 = _mm512_load_ps(&yim[0]);
                      _mm512_store_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_LT_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void clt_zmm16r4(const __m512 xre,
                                    const __m512 xim,
                                    const __m512 yre,
                                    const __m512 yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_ps_mask(xre,yre,_CMP_LT_OQ);
                         *eqi = _mm512_cmp_ps_mask(xim,yim,_CMP_LT_OQ);
              }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_zmm16r4_u(const float * __restrict xre,
                                      const float * __restrict xim,
                                      const float * __restrict yre,
                                      const float * __restrict yim,
                                      __mmask16 * __restrict eqr,
                                      __mmask16 * __restrict eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_loadu_ps(&xre[0]);
                      zmm1 = _mm512_loadu_ps(&yre[0]);
                      _mm512_storeu_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_NEQ_OQ));
                      zmm2 = _mm512_loadu_ps(&xim[0]);
                      zmm3 = _mm512_loadu_ps(&yim[0]);
                      _mm512_storeu_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_NEQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                      const float * __restrict __ATTR_ALIGN__(64) xim,
                                      const float * __restrict __ATTR_ALIGN__(64) yre,
                                      const float * __restrict __ATTR_ALIGN__(64) yim,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                      __mmask16 * __restrict __ATTR_ALIGN__(64) eqi ) {

                      register __m512 zmm0,zmm1,zmm2,zmm3;
                      zmm0 = _mm512_load_ps(&xre[0]);
                      zmm1 = _mm512_load_ps(&yre[0]);
                      _mm512_store_ps(&eqr[0],
                                       _mm512_cmp_ps_mask(zmm0,zmm1,_CMP_NEQ_OQ));
                      zmm2 = _mm512_load_ps(&xim[0]);
                      zmm3 = _mm512_load_ps(&yim[0]);
                      _mm512_store_ps(&eqi[0],
                                       _mm512_cmp_ps_mask(zmm2,zmm3,_CMP_NEQ_OQ));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cneq_zmm16r4(const __m512 xre,
                                    const __m512 xim,
                                    const __m512 yre,
                                    const __m512 yim,
                                    __mmask16 * __restrict eqr,
                                    __mmask16 * __restrict eqi) {

                         *eqr = _mm512_cmp_ps_mask(xre,yre,_CMP_NEQ_OQ);
                         *eqi = _mm512_cmp_ps_mask(xim,yim,_CMP_NEQ_OQ);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict cexpr,
                                       float * __restrict cexpi ) {

                        register const __m512 I = _mm512_set1_ps(1.0f);
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_loadu_ps(&xre[0]);
                        zmm1  = _mm512_loadu_ps(&xim[0]);
                        zmm2  = xexpf(zmm0);
                        zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                        _mm512_storeu_ps(&cexpr[0],zmm3);
                        zmm4  = _mm512_mul_ps(zmm2,_mm512_mul_ps(xsinf(zmm1),I));
                        _mm512_storeu_ps(&cexpi[0],zmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       float * __restrict __ATTR_ALIGN__(64) cexpr,
                                       float * __restrict __ATTR_ALIGN__(64) cexpi ) {

                        register const __m512 I = _mm512_set1_ps(1.0f);
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        zmm0  = _mm512_load_ps(&xre[0]);
                        zmm1  = _mm512_load_ps(&xim[0]);
                        zmm2  = xexpf(zmm0);
                        zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                        _mm512_store_ps(&cexpr[0],zmm3);
                        zmm4  = _mm512_mul_ps(zmm2,_mm512_mul_ps(xsinf(zmm1),I));
                        _mm512_store_ps(&cexpi[0],zmm4);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cexp_zmm16r4(const __m512 xre,
                                     const __m512 xim,
                                     __m512 * __restrict cexpr,
                                     __m512 * __restrict cexpi) {

                        register const __m512 I = _mm512_set1_ps(1.0f);
                        register __m512 zmm0;
                        zmm0   = xexpf(xre);
                        *cexpr = _mm512_mul_ps(zmm0,xcosf(xim));
                        *cexpi = _mm512_mul_ps(zmm0,_mm512_mul_ps(xsinf(xim),I));
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cexp_zmm16r4(const __m512 xre,
                                          const __m512 xim) {
                                     
                        zmm16c4_t cv;
                        register const __m512 I = _mm512_set1_ps(1.0f);
                        register __m512 zmm0;
                        zmm0   = xexpf(xre);
                        cv.re = _mm512_mul_ps(zmm0,xcosf(xim));
                        cv.im = _mm512_mul_ps(zmm0,_mm512_mul_ps(xsinf(xim),I));
                        return (cv);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_zmm16r4_u(const float * __restrict rho,
                                         const float * __restrict tht,
                                         float * __restrict  re,
                                         float * __restrict  im) {

                         register __m512 zmm0,zmm1,zmm2,zmm3;
                         zmm0 = _mm512_loadu_ps(&rho[0]);
                         zmm1 = _mm512_loadu_ps(&tht[0]);
                         zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                         _mm512_storeu_ps(&re[0],zmm2);
                         zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                         _mm512_storeu_ps(&im[0],zmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) rho,
                                         const float * __restrict __ATTR_ALIGN__(64) tht,
                                         float * __restrict  __ATTR_ALIGN__(64) re,
                                         float * __restrict  __ATTR_ALIGN__(64) im) {

                         register __m512 zmm0,zmm1,zmm2,zmm3;
                         zmm0 = _mm512_load_ps(&rho[0]);
                         zmm1 = _mm512_load_ps(&tht[0]);
                         zmm2 = _mm512_mul_ps(zmm0,xcosf(zmm1)); //tht
                         _mm512_store_ps(&re[0],zmm2);
                         zmm3 = _mm512_mul_ps(zmm0,xsinf(zmm1)); //tht
                         _mm512_store_ps(&im[0],zmm3);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cpolar_zmm16r4(const __m512 rho,
                                       const __m512 tht,
                                       __m512 * __restrict re,
                                       __m512 * __restrict im) {

                        register __m512 zmm0,zmm1;
                        zmm0 = _mm512_mul_ps(rho,xcosf(tht));
                        *re  = zmm0;
                        zmm1 = _mm512_mul_ps(rho,xsinf(tht));
                        *im  = zmm1;
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t cpolar_zmm16r4(const __m512 rho,
                                            const __m512 tht) {
                                      
                        zmm16c4_t cv
                        register __m512 zmm0,zmm1;
                        zmm0 = _mm512_mul_ps(rho,xcosf(tht));
                        cv.re  = zmm0;
                        zmm1 = _mm512_mul_ps(rho,xsinf(tht));
                        cv.im  = zmm1;
                        return (cv);
              }
              


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_zmm16r4_u(const float * __restrict xre,
                                       const float * __restrict xim,
                                       float * __restrict wrkc,
                                       float * __restrict csqr,
                                       float * __restrict csqi) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cabs_zmm16r4_u(xre,xim,wrkc);
                        zmm0  = _mm512_loadu_ps(&xre[0]);
                        zmm1  = _mm512_loadu_ps(&wrkc[0]);
                        zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                        _mm512_storeu_ps(&csqr[0],_mm512_sqrt_ps(zmm2));
                        zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                        _mm512_storeu_ps(&csqi[0],_mm512_sqrt_ps(zmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                       const float * __restrict __ATTR_ALIGN__(64) xim,
                                       float * __restrict __ATTR_ALIGN__(64) wrkc,
                                       float * __restrict __ATTR_ALIGN__(64) csqr,
                                       float * __restrict __ATTR_ALIGN__(64) csqi) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cabs_zmm16r4_a(xre,xim,wrkc);
                        zmm0  = _mm512_load_ps(&xre[0]);
                        zmm1  = _mm512_load_ps(&wrkc[0]);
                        zmm2  = _mm512_mul_ps(half,_mm512_add_ps(zmm1,zmm0));
                        _mm512_store_ps(&csqr[0],_mm512_sqrt_ps(zmm2));
                        zmm3  = _mm512_mul_ps(half,_mm512_sub_ps(zmm1,zmm0));
                        _mm512_store_ps(&csqi[0],_mm512_sqrt_ps(zmm3));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void csqrt_zmm16r4(const __m512 xre,
                                      const __m512 xim,
                                      __m512 * __restrict wrkc,
                                      __m512 * __restrict csqr,
                                      __m512 * __restrict csqi) {

                       register __m512 zmm0,zmm1;
                       const register __m512 half = _mm512_set1_ps(0.5f); 
                       cabs_zmm16r4(xre,xim,wrkc);
                       zmm0  = _mm512_mul_ps(half,_mm512_add_ps(*wrkc,xre));
                       *csqr = zmm0;
                       zmm1  = _mm512_mul_ps(half,_mm512_sub_ps(*wrkc,xre));
                       *csqi = zmm1; 
              }
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   zmm16c4_t csqrt_zmm16r4(const __m512 xre,
                                           const __m512 xim,
                                          __m512 * __restrict wrkc) {
                                          
                       zmm16c4_t cv;
                       register __m512 zmm0,zmm1;
                       const register __m512 half = _mm512_set1_ps(0.5f); 
                       cabs_zmm16r4(xre,xim,wrkc);
                       zmm0  = _mm512_mul_ps(half,_mm512_add_ps(*wrkc,xre));
                       cv.re = zmm0;
                       zmm1  = _mm512_mul_ps(half,_mm512_sub_ps(*wrkc,xre));
                       cv.im = zmm1; 
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_zmm16r4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict zre,
                                             float * __restrict zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        rep  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        zmm0 = _mm512_mul_ps(rep,rep);
                        zmm1 = _mm512_mul_ps(imp,imp);
                        zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                        _mm512_storeu_ps(&zre[0], _mm512_div_ps(rep,zmm2));
                        _mm512_storeu_ps(&zim[0], _mm512_div_ps(imp,zmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                             const float * __restrict __ATTR_ALIGN__(64) xim,
                                             const float * __restrict __ATTR_ALIGN__(64) yre,
                                             const float * __restrict __ATTR_ALIGN__(64) yim,
                                             float * __restrict __ATTR_ALIGN__(64) zre,
                                             float * __restrict __ATTR_ALIGN__(64) zim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        rep  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        imp  = _mm512_fmadd_pd(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        zmm0 = _mm512_mul_ps(rep,rep);
                        zmm1 = _mm512_mul_ps(imp,imp);
                        zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                        _mm512_store_ps(&zre[0], _mm512_div_ps(rep,zmm2));
                        _mm512_store_ps(&zim[0], _mm512_div_ps(imp,zmm2));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_prod_zmm16r4(  const __m512  xre,
                                             const __m512  xim,
                                             const __m512  yre,
                                             const __m512  yim,
                                             __m512 * __restrict zre,
                                             __m512 * __restrict zim) {

                        register __m512 rep,imp,zmm0,zmm1,zmm2;
                        rep  = _mm512_fmsub_ps(xre,yre,
                                               _mm512_mul_ps(xim,yim));
                        imp  = _mm512_fmadd_pd(xim,yre,
                                               _mm512_mul_ps(xre,yim));
                        zmm0 = _mm512_mul_ps(rep,rep);
                        zmm1 = _mm512_mul_ps(imp,imp);
                        zmm2 = _mm512_sqrt_ps(_mm512_add_ps(zmm0,zmm1));
                        *zre = _mm512_div_ps(rep,zmm2);
                        *zim = _mm512_div_ps(imp,zmm2);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_zmm16r4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        sre  = _mm512_reduce_add_ps(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        sim  = _mm512_reduce_add_ps(imp);
                        *mim = sim*inv16;
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                             const float * __restrict __ATTR_ALIGN__(64) xim,
                                             const float * __restrict __ATTR_ALIGN__(64) yre,
                                             const float * __restrict __ATTR_ALIGN__(64) yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        sre = 0.0f;
                        rep  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        sre  = _mm512_reduce_add_ps(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        sim  = _mm512_reduce_add_ps(imp);
                        *mim = sim*inv16;
              } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_prod_zmm16r4(const __m512 xre,
                                           const __m512 xim,
                                           const __m512 yre,
                                           const __m512 yim,
                                           float * __restrict mre,
                                           float * __restrict min) {

                        register __m512 rep,imp;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre = 0.0f;
                        rep  = _mm512_fmsub_ps(xre,yre,
                                               _mm512_mul_ps(xim,yim));
                        sre  = _mm512_reduce_add_ps(rep);
                        *mre = sre*inv16;
                        sim  = 0.0f;
                        imp  = _mm512_fmadd_ps(xim,yre,
                                               _mm512_mul_ps(xre,yim));
                        sim  = _mm512_reduce_add_ps(imp);
                        *mim = sim*inv16;
             }


               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_zmm16r4_u(const float * __restrict xre,
                                             const float * __restrict xim,
                                             const float * __restrict yre,
                                             const float * __restrict yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        imp  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        sim  = 0.0f;
                        den  = _mm512_fmadd_ps(zmm1,zmm1,
                                               _mm512_mul_ps(zmm3,zmm3));
                        rquot = _mm512_div_ps(rep,den);
                        sre   = _mm512_reduce_add_ps(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm512_div_ps(imp,den);
                        sim   = _mm512_reduce_add_ps(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                             const float * __restrict __ATTR_ALIGN__(64) xim,
                                             const float * __restrict __ATTR_ALIGN__(64) yre,
                                             const float * __restrict __ATTR_ALIGN__(64) yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        sre  = 0.0f;
                        rep  = _mm512_fmsub_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        imp  = _mm512_fmadd_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        sim  = 0.0f;
                        den  = _mm512_fmadd_ps(zmm1,zmm1,
                                               _mm512_mul_ps(zmm3,zmm3));
                        rquot = _mm512_div_ps(rep,den);
                        sre   = _mm512_reduce_add_ps(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm512_div_ps(imp,den);
                        sim   = _mm512_reduce_add_ps(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_quot_zmm16r4(  const __m512 xre,
                                             const __m512 xim,
                                             const __m512 yre,
                                             const __m512 yim,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m512 rep,imp,den,rquot,iquot;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre  = 0.0f;
                        rep  = _mm512_fmsub_ps(xre,yre,
                                               _mm512_mul_ps(xim,yim));
                        imp  = _mm512_fmadd_ps(xim,yre,
                                               _mm512_mul_ps(xre,yim));
                        sim  = 0.0f;
                        den  = _mm512_fmadd_ps(yre,yre,
                                               _mm512_mul_ps(yim,yim));
                        rquot = _mm512_div_ps(rep,den);
                        sre   = _mm512_reduce_add_ps(rquot);
                        *mre  = sre*inv16;
                        iquot = _mm512_div_ps(imp,den);
                        sim   = _mm512_reduce_add_ps(iquot);
                        *mim  = sre*inv16;
              }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_zmm16r4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim ) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp,magc1,magc2,vcmag;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        rep  = _mm512_fmadd_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        magc1= _mm512_mul_ps(rep,rep);
                        imp  = _mm512_fmsub_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        magc2= _mm512_mul_ps(imp,imp);
                        vcmag= _mm512_sqrt_ps(_mm512_add_ps(magc1,magc2));
                        _mm512_storeu_ps(&mre[0], _mm512_div_ps(rep,vcmag));
                        _mm512_storeu_ps(&mim[0], _mm512_div_ps(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                              const float * __restrict __ATTR_ALIGN__(64) xim,
                                              const float * __restrict __ATTR_ALIGN__(64) yre,
                                              const float * __restrict __ATTR_ALIGN__(64) yim,
                                              float * __restrict __ATTR_ALIGN__(64) mre,
                                              float * __restrict __ATTR_ALIGN__(64) mim ) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 rep,imp,magc1,magc2,vcmag;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        rep  = _mm512_fmad_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        magc1= _mm512_mul_ps(rep,rep);
                        imp  = _mm512_fmsub_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        magc2= _mm512_mul_ps(imp,imp);
                        vcmag= _mm512_sqrt_ps(_mm512_add_ps(magc1,magc2));
                        _mm512_store_ps(&mre[0], _mm512_div_ps(rep,vcmag));
                        _mm512_store_ps(&mim[0], _mm512_div_ps(imp,vcmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnorm_cprod_zmm16r4(const __m512 xre,
                                            const __m512 xim,
                                            const __m512 yre,
                                            const __m512 yim,
                                            __m512 * __restrict mre,
                                            __m512 * __restrict mim) {

                        register __m512 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm512_fmad_ps(xre,yre,
                                               _mm512_mul_ps(xim,yim));
                        magc1= _mm512_mul_ps(rep,rep);
                        imp  = _mm512_fmsub_ps(xim,yre,
                                               _mm512_mul_ps(xre,yim));
                        magc2= _mm512_mul_ps(imp,imp);
                        vcmag= _mm512_sqrt_ps(_mm512_add_ps(magc1,magc2));
                        *mre = _mm512_div_ps(rep,vcmag);
                        *mim = _mm512_div_ps(imp,vcmag)
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_zmm16r4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim) {
                      
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        re   = _mm512_fmadd_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        sre  = _mm512_reduce_add_ps(re);
                        *mre = sre*inv16;
                        im   = _mm512_fmsub_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        sim  = _mm512_reduce_add_ps(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                              const float * __restrict __ATTR_ALIGN__(64) xim,
                                              const float * __restrict __ATTR_ALIGN__(64) yre,
                                              const float * __restrict __ATTR_ALIGN__(64) yim,
                                              float * __restrict mre,
                                              float * __restrict mim) {
                      
                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        re   = _mm512_fmadd_ps(zmm0,zmm1,
                                               _mm512_mul_ps(zmm2,zmm3));
                        sre  = _mm512_reduce_add_ps(re);
                        *mre = sre*inv16;
                        im   = _mm512_fmsub_ps(zmm2,zmm1,
                                               _mm512_mul_ps(zmm0,zmm3));
                        sim  = _mm512_reduce_add_ps(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmean_cprod_zmm16r4(const __m512 xre,
                                            const __m512 xim,
                                            const __m512 yre,
                                            const __m512 yim,
                                            float * __restrict mre,
                                            float * __restrict min) {

                        register __m512 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm512_fmadd_ps(xre,yre,
                                               _mm512_mul_ps(xim,yim));
                        sre  = _mm512_reduce_add_ps(re);
                        *mre = sre*inv16;
                        im   = _mm512_fmsub_ps(xim,yre,
                                               _mm512_mul_ps(xre,yim));
                        sim  = _mm512_reduce_add_ps(im);
                        *mim = sim*inv16;
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_zmm16r4_u(const float * __restrict xre,
                                              const float * __restrict xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        register __m512 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm512_loadu_ps(&xre[0]);
                        sre  = _mm512_reduce_add_ps(re);
                        *mre = sre*inv16;
                        im   = _mm512_loadu_ps(&xim[0]);
                        sim  = _mm512_reduce_add_ps(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                              const float * __restrict __ATTR_ALIGN__(64) xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        register __m512 re,im;
                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        re   = _mm512_load_ps(&xre[0]);
                        sre  = _mm512_reduce_add_ps(re);
                        *mre = sre*inv16;
                        im   = _mm512_load_ps(&xim[0]);
                        sim  = _mm512_reduce_add_ps(im);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void arith_cmean_zmm16r4(  const __m512 xre,
                                              const __m512 xim,
                                              float * __restrict mre,
                                              float * __restrict min) {

                        constexpr float inv16 = 0.0625f;
                        float sre,sim;
                        sre  = _mm512_reduce_add_ps(xre);
                        *mre = sre*inv16;
                        sim  = _mm512_reduce_add_ps(xim);
                        *mim = sim*inv16; 
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_zmm16r4_u( const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre,
                                              float * __restrict mim ) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 re,im,cvmag;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        cvmag= _mm512_sqrt_ps(_mm512_fmadd_ps(zmm0,zmm1,
                                                              _mm512_mul_ps(zmm2,zmm3)));
                        _mm512_storeu_ps(&mre[0], _mm512_div_ps(zmm0,cvmag));
                        _mm512_storeu_ps(&mim[0], _mm512_div_ps(zmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                              const float * __restrict __ATTR_ALIGN__(64) xim,
                                              const float * __restrict __ATTR_ALIGN__(64) yre,
                                              const float * __restrict __ATTR_ALIGN__(64) yim,
                                              float * __restrict __ATTR_ALIGN__(64) mre,
                                              float * __restrict __ATTR_ALIGN__(64) mim ) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 re,im,cvmag;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        cvmag= _mm512_sqrt_ps(_mm512_fmadd_ps(zmm0,zmm1,
                                                              _mm512_mul_ps(zmm2,zmm3)));
                        _mm512_store_ps(&mre[0], _mm512_div_ps(zmm0,cvmag));
                        _mm512_store_ps(&mim[0], _mm512_div_ps(zmm2,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cnormalize_zmm16r4( const __m512 xre,
                                            const __m512 xim,
                                            const __m512 yre,
                                            const __m512 yim,
                                            __m512 * __restrict mre,
                                            __m512 * __restrict mim ) {

                        register __m512 re,im,cvmag;
                        cvmag= _mm512_sqrt_ps(_mm512_fmadd_ps(xre,yre,
                                                    _mm512_mul_ps(xim,yim)));
                        *mre = _mm512_div_ps(xre,cvmag));
                        *mim =  _mm512_div_ps(xim,cvmag));
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_zmm16r4_u( const float * __restrict xre,
                                              const float * __restrict xim,
                                              const float * __restrict yre,
                                              const float * __restrict yim,
                                              float * __restrict mre) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 cvmag;
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&yre[0]);
                        zmm2 = _mm512_loadu_ps(&xim[0]);
                        zmm3 = _mm512_loadu_ps(&yim[0]);
                        cvmag= _mm512_sqrt_ps(_mm512_fmadd_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3)));
                        _mm512_storeu_ps(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                              const float * __restrict __ATTR_ALIGN__(64) xim,
                                              const float * __restrict __ATTR_ALIGN__(64) yre,
                                              const float * __restrict __ATTR_ALIGN__(64) yim,
                                              float * __restrict __ATTR_ALIGN__(64) mre) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 cvmag;
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&yre[0]);
                        zmm2 = _mm512_load_ps(&xim[0]);
                        zmm3 = _mm512_load_ps(&yim[0]);
                        cvmag= _mm512_sqrt_ps(_mm512_fmadd_ps(zmm0,zmm1,
                                                          _mm512_mul_ps(zmm2,zmm3)));
                        _mm512_store_ps(&mre[0], cvmag);
             }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cmagnitude_zmm16r4(   const __m512 xre,
                                              const __m512 xim,
                                              const __m512 yre,
                                              const __m512 yim,
                                              __m512 * __restrict  mre) {

                        register __m512 cvmag;
                        cvmag= _mm512_sqrt_ps(_mm512_fmadd_ps(xre,yre,
                                                          _mm512_mul_ps(xim,yim)));
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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_2xr8_c8_unroll16x(double * __restrict __ATTR_ALIGN__(64) xre,
                                               double * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(64) vc,
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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
                   static inline
                   void copy_c8_2xr8_unroll16x( double * __restrict __ATTR_ALIGN__(64) xre,
                                                double * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<double> * __restrict __ATTR_ALIGN__(64) vc,
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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_zmm16r4_u(const float * __restrict xre,
                                               const float * __restrict xim,
                                               std::complex<float> * __restrict yc) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 zmm4,zmm5;
                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        zmm0 = _mm512_loadu_ps(&xre[0]);
                        zmm1 = _mm512_loadu_ps(&xim[0]);
                        zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                        zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                        zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                        _mm512_storeu_ps(&pyc[0],zmm4);
                        zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                        _mm512_storeu_ps(&pyc[16],zmm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_zmm16r4(  const __m512 xre,
                                               const __m512 xim,
                                               std::complex<float> * __restrict yc) {

                        register __m512 zmm2,zmm3;
                        register __m512 zmm4,zmm5;
                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        zmm2 = _mm512_unpacklo_ps(xre,xim);
                        zmm3 = _mm512_unpackhi_ps(xre,xim);
                        zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                        _mm512_storeu_ps(&pyc[0],zmm4);
                        zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                        _mm512_storeu_ps(&pyc[16],zmm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void copy_2xr4_c4_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                               const float * __restrict __ATTR_ALIGN__(64) xim,
                                               std::complex<float> * __restrict __ATTR_ALIGN__(64) yc) {

                        register __m512 zmm0,zmm1,zmm2,zmm3;
                        register __m512 zmm4,zmm5;
                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        zmm0 = _mm512_load_ps(&xre[0]);
                        zmm1 = _mm512_load_ps(&xim[0]);
                        zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                        zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                        zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                        _mm512_store_ps(&pyc[0],zmm4);
                        zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                        _mm512_store_ps(&pyc[16],zmm5);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
                   void copy_2r4c4_zmm16r4_blocked_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   std::complex<float> * __restrict __ATTR_ALIGN__(64) yc,
                                                   const int32_t n) { // size of array std::complex<float> len elements.

                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        if(n <= 1) {
                            return;
                        }
                        else if(n <= 16) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_load_ps(&xre[0]);
                            zmm1 = _mm512_load_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[16],zmm5);
                        }
                        else if(n <= 32) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_load_ps(&xre[0]);
                            zmm1 = _mm512_load_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[16]);
                            zmm1 = _mm512_load_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[48],zmm5);
                       }
                       else if(n <= 64) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_load_ps(&xre[0]);
                            zmm1 = _mm512_load_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[16]);
                            zmm1 = _mm512_load_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[32]);
                            zmm1 = _mm512_load_ps(&xim[32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[48]);
                            zmm1 = _mm512_load_ps(&xim[48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[112],zmm5);
                       }
                       else if(n <= 128) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_load_ps(&xre[0]);
                            zmm1 = _mm512_load_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[16]);
                            zmm1 = _mm512_load_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[32]);
                            zmm1 = _mm512_load_ps(&xim[32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[48]);
                            zmm1 = _mm512_load_ps(&xim[48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[112],zmm5);
                            zmm0 = _mm512_load_ps(&xre[64]);
                            zmm1 = _mm512_load_ps(&xim[64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[144],zmm5);
                            zmm0 = _mm512_load_ps(&xre[80]);
                            zmm1 = _mm512_load_ps(&xim[80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[176],zmm5);
                            zmm0 = _mm512_load_ps(&xre[96]);
                            zmm1 = _mm512_load_ps(&xim[96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[208],zmm5);
                            zmm0 = _mm512_load_ps(&xre[112]);
                            zmm1 = _mm512_load_ps(&xim[112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[240],zmm5);
                       }
                       else if(n <= 256) {
#if (GMS_MAN_PREFETCH) == 1
                               _mm_prefetch((const char*)&xre[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[128],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[128],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[160],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[160],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[192],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[192],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[224],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[224],_MM_HINT_T0);
#endif                           
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_load_ps(&xre[0]);
                            zmm1 = _mm512_load_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[16]);
                            zmm1 = _mm512_load_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[32]);
                            zmm1 = _mm512_load_ps(&xim[32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[48]);
                            zmm1 = _mm512_load_ps(&xim[48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[112],zmm5);
                            zmm0 = _mm512_load_ps(&xre[64]);
                            zmm1 = _mm512_load_ps(&xim[64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[144],zmm5);
                            zmm0 = _mm512_load_ps(&xre[80]);
                            zmm1 = _mm512_load_ps(&xim[80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[176],zmm5);
                            zmm0 = _mm512_load_ps(&xre[96]);
                            zmm1 = _mm512_load_ps(&xim[96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[208],zmm5);
                            zmm0 = _mm512_load_ps(&xre[112]);
                            zmm1 = _mm512_load_ps(&xim[112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[240],zmm5);
                            zmm0 = _mm512_load_ps(&xre[128]);
                            zmm1 = _mm512_load_ps(&xim[128]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[256],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[272],zmm5);
                            zmm0 = _mm512_load_ps(&xre[144]);
                            zmm1 = _mm512_load_ps(&xim[144]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[288],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[304],zmm5);
                            zmm0 = _mm512_load_ps(&xre[160]);
                            zmm1 = _mm512_load_ps(&xim[160]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[320],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[336],zmm5);
                            zmm0 = _mm512_load_ps(&xre[176]);
                            zmm1 = _mm512_load_ps(&xim[176]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[352],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[368],zmm5);
                            zmm0 = _mm512_load_ps(&xre[192]);
                            zmm1 = _mm512_load_ps(&xim[192]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[384],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[400],zmm5);
                            zmm0 = _mm512_load_ps(&xre[208]);
                            zmm1 = _mm512_load_ps(&xim[208]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[416],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[432],zmm5);
                            zmm0 = _mm512_load_ps(&xre[224]);
                            zmm1 = _mm512_load_ps(&xim[224]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[448],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[464],zmm5);
                            zmm0 = _mm512_load_ps(&xre[240]);
                            zmm1 = _mm512_load_ps(&xim[240]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[480],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[496],zmm5);
                       }
                      else if(n > 256) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            int32_t i;
                       for(i = 0; (i+255) < n; i += 256) {
#if (GMS_MAN_PREFETCH) == 1
                            _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+224],_MM_HINT_T0);
#endif                         
                            zmm0 = _mm512_load_ps(&xre[i+0]);
                            zmm1 = _mm512_load_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+16]);
                            zmm1 = _mm512_load_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+32]);
                            zmm1 = _mm512_load_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+48]);
                            zmm1 = _mm512_load_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+112],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+64]);
                            zmm1 = _mm512_load_ps(&xim[i+64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+144],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+80]);
                            zmm1 = _mm512_load_ps(&xim[i+80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+176],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+96]);
                            zmm1 = _mm512_load_ps(&xim[i+96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+208],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+112]);
                            zmm1 = _mm512_load_ps(&xim[i+112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+240],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+128]);
                            zmm1 = _mm512_load_ps(&xim[i+128]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+256],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+272],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+144]);
                            zmm1 = _mm512_load_ps(&xim[i+144]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+288],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+304],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+160]);
                            zmm1 = _mm512_load_ps(&xim[i+160]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+320],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+336],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+176]);
                            zmm1 = _mm512_load_ps(&xim[i+176]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+352],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+368],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+192]);
                            zmm1 = _mm512_load_ps(&xim[i+192]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+384],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+400],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+208]);
                            zmm1 = _mm512_load_ps(&xim[i+208]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+416],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+432],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+224]);
                            zmm1 = _mm512_load_ps(&xim[i+224]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+448],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+464],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+240]);
                            zmm1 = _mm512_load_ps(&xim[i+240]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+480],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+496],zmm5);
                     }

                       for(; (i+191) < n; i += 192) {
                            zmm0 = _mm512_load_ps(&xre[i+0]);
                            zmm1 = _mm512_load_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+16]);
                            zmm1 = _mm512_load_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+32]);
                            zmm1 = _mm512_load_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+48]);
                            zmm1 = _mm512_load_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+112],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+64]);
                            zmm1 = _mm512_load_ps(&xim[i+64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+144],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+80]);
                            zmm1 = _mm512_load_ps(&xim[i+80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+176],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+96]);
                            zmm1 = _mm512_load_ps(&xim[i+96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+208],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+112]);
                            zmm1 = _mm512_load_ps(&xim[i+112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+240],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+128]);
                            zmm1 = _mm512_load_ps(&xim[i+128]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+256],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+272],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+144]);
                            zmm1 = _mm512_load_ps(&xim[i+144]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+288],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+304],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+160]);
                            zmm1 = _mm512_load_ps(&xim[i+160]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+320],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+336],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+176]);
                            zmm1 = _mm512_load_ps(&xim[i+176]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+352],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+368],zmm5);
                     }    

                       for(; (i+127) < n; i += 128) {
                            zmm0 = _mm512_load_ps(&xre[i+0]);
                            zmm1 = _mm512_load_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+16]);
                            zmm1 = _mm512_load_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+32]);
                            zmm1 = _mm512_load_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+48]);
                            zmm1 = _mm512_load_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+112],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+64]);
                            zmm1 = _mm512_load_ps(&xim[i+64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+144],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+80]);
                            zmm1 = _mm512_load_ps(&xim[i+80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+176],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+96]);
                            zmm1 = _mm512_load_ps(&xim[i+96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+208],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+112]);
                            zmm1 = _mm512_load_ps(&xim[i+112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+240],zmm5); 
                     }

                       for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_load_ps(&xre[i+0]);
                            zmm1 = _mm512_load_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+16]);
                            zmm1 = _mm512_load_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+32]);
                            zmm1 = _mm512_load_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+48]);
                            zmm1 = _mm512_load_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+112],zmm5);
                     }

                       for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_load_ps(&xre[i+0]);
                            zmm1 = _mm512_load_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_load_ps(&xre[i+16]);
                            zmm1 = _mm512_load_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+48],zmm5);
                     }

                       for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_load_ps(&xre[i+0]);
                            zmm1 = _mm512_load_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+16],zmm5);
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
                   void copy_2r4c4_zmm16r4_blocked_u(const float * __restrict  xre,
                                                   const float * __restrict  xim,
                                                   std::complex<float> * __restrict  yc,
                                                   const int32_t n) { // size of array std::complex<float> len elements.

                        float * __restrict pyc = reinterpret_cast<float*>(&yc[0]);
                        if(n <= 1) {
                            return;
                        }
                        else if(n <= 16) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_loau_ps(&xre[0]);
                            zmm1 = _mm512_loadu_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[16],zmm5);
                        }
                        else if(n <= 32) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_loadu_ps(&xre[0]);
                            zmm1 = _mm512_loadu_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[16]);
                            zmm1 = _mm512_loadu_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[48],zmm5);
                       }
                       else if(n <= 64) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_loadu_ps(&xre[0]);
                            zmm1 = _mm512_loadu_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[16]);
                            zmm1 = _mm512_loadu_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[32]);
                            zmm1 = _mm512_loadu_ps(&xim[32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[48]);
                            zmm1 = _mm512_loadu_ps(&xim[48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[112],zmm5);
                       }
                       else if(n <= 128) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_loadu_ps(&xre[0]);
                            zmm1 = _mm512_loadu_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[16]);
                            zmm1 = _mm512_loadu_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[32]);
                            zmm1 = _mm512_loadu_ps(&xim[32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[48]);
                            zmm1 = _mm512_loadu_ps(&xim[48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[112],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[64]);
                            zmm1 = _mm512_loadu_ps(&xim[64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[144],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[80]);
                            zmm1 = _mm512_loadu_ps(&xim[80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[176],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[96]);
                            zmm1 = _mm512_loadu_ps(&xim[96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[208],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[112]);
                            zmm1 = _mm512_loadu_ps(&xim[112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[240],zmm5);
                       }
                       else if(n <= 256) {
#if (GMS_MAN_PREFETCH) == 1
                               _mm_prefetch((const char*)&xre[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[0],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[32],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[64],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[96],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[128],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[128],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[160],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[160],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[192],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[192],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xre[224],_MM_HINT_T0);
                               _mm_prefetch((const char*)&xim[224],_MM_HINT_T0);
#endif                           
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            zmm0 = _mm512_loadu_ps(&xre[0]);
                            zmm1 = _mm512_loadu_ps(&xim[0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[16]);
                            zmm1 = _mm512_loadu_ps(&xim[16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[32]);
                            zmm1 = _mm512_loadu_ps(&xim[32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[48]);
                            zmm1 = _mm512_loadu_ps(&xim[48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[112],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[64]);
                            zmm1 = _mm512_loadu_ps(&xim[64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[144],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[80]);
                            zmm1 = _mm512_loadu_ps(&xim[80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[176],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[96]);
                            zmm1 = _mm512_loadu_ps(&xim[96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[208],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[112]);
                            zmm1 = _mm512_loadu_ps(&xim[112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[240],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[128]);
                            zmm1 = _mm512_loadu_ps(&xim[128]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[256],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[272],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[144]);
                            zmm1 = _mm512_loadu_ps(&xim[144]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[288],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[304],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[160]);
                            zmm1 = _mm512_loadu_ps(&xim[160]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[320],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[336],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[176]);
                            zmm1 = _mm512_loadu_ps(&xim[176]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[352],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[368],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[192]);
                            zmm1 = _mm512_loadu_ps(&xim[192]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[384],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[400],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[208]);
                            zmm1 = _mm512_loadu_ps(&xim[208]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[416],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[432],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[224]);
                            zmm1 = _mm512_loadu_ps(&xim[224]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[448],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[464],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[240]);
                            zmm1 = _mm512_loadu_ps(&xim[240]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[480],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[496],zmm5);
                       }
                      else if(n > 256) {
                            register __m512 zmm0,zmm1,zmm2,zmm3;
                            register __m512 zmm4,zmm5;
                            int32_t i;
                       for(i = 0; (i+255) < n; i += 256) {
#if (GMS_MAN_PREFETCH) == 1
                            _mm_prefetch((const char*)&xre[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+0],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+192],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xre[i+224],_MM_HINT_T0);
                            _mm_prefetch((const char*)&xim[i+224],_MM_HINT_T0);
#endif                         
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+16]);
                            zmm1 = _mm512_loadu_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+32]);
                            zmm1 = _mm512_loadu_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+48]);
                            zmm1 = _mm512_loadu_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+112],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+64]);
                            zmm1 = _mm512_loadu_ps(&xim[i+64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+144],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+80]);
                            zmm1 = _mm512_loadu_ps(&xim[i+80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+176],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+96]);
                            zmm1 = _mm512_loadu_ps(&xim[i+96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+208],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+112]);
                            zmm1 = _mm512_loadu_ps(&xim[i+112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+240],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+128]);
                            zmm1 = _mm512_loadu_ps(&xim[i+128]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+256],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+272],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+144]);
                            zmm1 = _mm512_loadu_ps(&xim[i+144]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+288],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+304],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+160]);
                            zmm1 = _mm512_loadu_ps(&xim[i+160]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+320],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+336],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+176]);
                            zmm1 = _mm512_loadu_ps(&xim[i+176]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+352],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+368],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+192]);
                            zmm1 = _mm512_loadu_ps(&xim[i+192]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+384],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+400],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+208]);
                            zmm1 = _mm512_loadu_ps(&xim[i+208]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+416],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+432],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+224]);
                            zmm1 = _mm512_loadu_ps(&xim[i+224]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+448],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+464],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+240]);
                            zmm1 = _mm512_loadu_ps(&xim[i+240]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+480],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+496],zmm5);
                     }

                       for(; (i+191) < n; i += 192) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+16]);
                            zmm1 = _mm512_loadu_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+32]);
                            zmm1 = _mm512_loadu_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+48]);
                            zmm1 = _mm512_loadu_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+112],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+64]);
                            zmm1 = _mm512_loadu_ps(&xim[i+64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+144],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+80]);
                            zmm1 = _mm512_loadu_ps(&xim[i+80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+176],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+96]);
                            zmm1 = _mm512_loadu_ps(&xim[i+96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+208],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+112]);
                            zmm1 = _mm512_loadu_ps(&xim[i+112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+240],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+128]);
                            zmm1 = _mm512_loadu_ps(&xim[i+128]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+256],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+272],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+144]);
                            zmm1 = _mm512_loadu_ps(&xim[i+144]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+288],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+304],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+160]);
                            zmm1 = _mm512_loadu_ps(&xim[i+160]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+320],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+336],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+176]);
                            zmm1 = _mm512_loadu_ps(&xim[i+176]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+352],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+368],zmm5);
                     }    

                       for(; (i+127) < n; i += 128) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+16]);
                            zmm1 = _mm512_loadu_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+32]);
                            zmm1 = _mm512_loadu_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+48]);
                            zmm1 = _mm512_loadu_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+112],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+64]);
                            zmm1 = _mm512_loadu_ps(&xim[i+64]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+128],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+144],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+80]);
                            zmm1 = _mm512_loadu_ps(&xim[i+80]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+160],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+176],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+96]);
                            zmm1 = _mm512_loadu_ps(&xim[i+96]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_store_ps(&pyc[i+192],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_store_ps(&pyc[i+208],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+112]);
                            zmm1 = _mm512_loadu_ps(&xim[i+112]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+224],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+240],zmm5); 
                     }

                       for(; (i+63) < n; i += 64) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+16]);
                            zmm1 = _mm512_loadu_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+48],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+32]);
                            zmm1 = _mm512_loadu_ps(&xim[i+32]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+64],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+80],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+48]);
                            zmm1 = _mm512_loadu_ps(&xim[i+48]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+96],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+112],zmm5);
                     }

                       for(; (i+31) < n; i += 32) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+16],zmm5);
                            zmm0 = _mm512_loadu_ps(&xre[i+16]);
                            zmm1 = _mm512_loadu_ps(&xim[i+16]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+32],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+48],zmm5);
                     }

                       for(; (i+15) < n; i += 16) {
                            zmm0 = _mm512_loadu_ps(&xre[i+0]);
                            zmm1 = _mm512_loadu_ps(&xim[i+0]);
                            zmm2 = _mm512_unpacklo_ps(zmm0,zmm1);
                            zmm3 = _mm512_unpackhi_ps(zmm0,zmm1);
                            zmm4 = _mm512_shuffle_f32x4(zmm2,zmm3,0x0);
                            _mm512_storeu_ps(&pyc[i+0],zmm4);
                            zmm5 = _mm512_shuffle_f32x4(zmm2,zmm3,0x3);
                            _mm512_storeu_ps(&pyc[i+16],zmm5);
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
