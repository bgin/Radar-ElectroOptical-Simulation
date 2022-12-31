
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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
                   void cadd_zmm8r8(const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd yre,
                                     const __m512ddd yim,
                                     __m512ddd &     zre,
                                     __m512ddd &     zim) {
                     
                        register __m512ddd zmm0,zmm1;
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
                   void cadd_zmm8r8(const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd s,
                                     __m512ddd &     zre,
                                     __m512ddd &     zim) {

                        zre = _mm512_add_pd(xre,s);
                        zim = _mm512_add_pd(xim,s);
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
                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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
                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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
                   void csub_zmm8r8(const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd yre,
                                     const __m512ddd yim,
                                     __m512ddd &     zre,
                                     __m512ddd &     zim) {
                     
                        register __m512ddd zmm0,zmm1;
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
                   void csub_zmm8r8(const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd s,
                                     __m512ddd &     zre,
                                     __m512ddd &     zim) {

                        zre = _mm512_sub_pd(xre,s);
                        zim = _mm512_sub_pd(xim,s);
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
                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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
                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
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

                           register __m512ddd zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
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

                           register __m512ddd zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
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
                   void cmul_zmm8r8(const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd yre,
                                     const __m512ddd yim,
                                     __m512ddd &     zre,
                                     __m512ddd &     zim) {

                         register __m512ddd zmm0,zmm1;
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
                   void cmul_zmm8r8(const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd s,
                                     __m512ddd &     zre,
                                     __m512ddd &     zim) {

                        zre = _mm512_mul_pd(xre,s);
                        zim = _mm512_mul_pd(xim,s);
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

                           register __m512ddd zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
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

                           register __m512ddd zmm0,zmm1,zmm2,zmm3,zmm4,zmm5;
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3; 
                        register __m512ddd zmm4,zmm5,zmm6;
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3; 
                        register __m512ddd zmm4,zmm5,zmm6;
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
                   void cdiv_zmm8r8( const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd yre,
                                     const __m512ddd yim,
                                     __m512dd & zre,
                                     __m512dd & zim) {

                      register __m512ddd zmm0,zmm1,zmm2;
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
                   void cdiv_zmm8r8( const __m512ddd xre,
                                     const __m512ddd xim,
                                     const __m512ddd s,
                                     __m512ddd & zre,
                                     __m512ddd & zim) {

                        zre = _mm512_div_pd(xre,s);
                        zim = _mm512_div_pd(xim,s);
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3; 
                        register __m512ddd zmm4,zmm5,zmm6;
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

                        register __m512ddd zmm0,zmm1,zmm2,zmm3; 
                        register __m512ddd zmm4,zmm5,zmm6;
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
                   void cdiv_smith_zmm8r8_u( const double * __restrict xre,
                                             const double * __restrict xim,
                                             const double * __restrict yre,
                                             const double * __restrict yim,
                                             double *       __restrict zre,
                                             double *       __restrict zim) {

                        register __m512ddd zmm0,zmm1,zmm2,zmm3;
                        register __m512ddd r,den;
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

                        register __m512dd zmm0,zmm1,zmm2,zmm3;
                        register __m512dd r,den;
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
                   void cdiv_smith_zmm8r8( const __m512dd xre,
                                           const __m512dd xim,
                                           const __m512dd yre,
                                           const __m512dd yim,
                                           __m512dd * __restrict zre,
                                           __m512dd * __restrict zim) {

                        register __m512dd r,den;
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
                   void cconj_zmm8r8_u(const double * __restrict im,
                                        double * __restrict  conj) {

                        register __m512d zmm0;
                        const register __m512d none = _mm512_set1_pd(-1.0f);
                        zmm0 = _mm512_loadu_pd(&im[0]);
                        _mm512_storeu_pd(&conj[0], _mm512_mul_pd(none,zmm0));
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void cconj_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) im,
                                        double * __restrict  __ATTR_ALIGN__(64) conj) {

                        register __m512d zmm0;
                        const register __m512d none = _mm512_set1_pd(-1.0f);
                        zmm0 = _mm512_load_pd(&im[0]);
                        _mm512_store_pd(&conj[0], _mm512_mul_pd(none,zmm0));
               }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d cconj_zmm8r8(const __m512d im) {
                          
                         const register __m512d none = _mm512_set1_pd(-1.0f);
                         register __m512d conj;
                         conj = _mm512_mul_pd(none,im);
                         return (conj); 
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


                                         


      } // math


} // gms















#endif /*__GMS_COMPLEX_ZMM8R8_HPP__*/
