
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

                        register __m512d zmm0,zmm1,zmm2,zmm3;
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

                        register __m512d zmm0,zmm1,zmm2,zmm3;
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



                                         


      } // math


} // gms















#endif /*__GMS_COMPLEX_ZMM8R8_HPP__*/
