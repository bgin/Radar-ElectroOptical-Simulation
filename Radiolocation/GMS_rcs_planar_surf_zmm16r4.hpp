

#ifndef __GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__
#define __GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__ 180420231543



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

    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_FULLVER =
      1000U*GMS_RCS_PLANAR_SURF_ZMM16R4_MAJOR+
      100U*GMS_RCS_PLANAR_SURF_ZMM16R4_MINOR+
      10U*GMS_RCS_PLANAR_SURF_ZMM16R4_MICRO;
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_CREATION_DATE = "18-04-2023 15:43 PM +00200 (TUE 18 APR 2023 GMT+2)";
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_DESCRIPTION   = "AVX512 optimized Planar Surfaces Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"



namespace  gms {


         namespace radiolocation {


                 /*
                       Complex impedances.
                       Formula 7.1-6
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm16r4(const __m512 tht,
                                        const __m512 mur,
                                        const __m512 mui,
                                        const __m512 epsr,
                                        const __m512 epsi,
                                        __m512 * __restrict zr,
                                        __m512 * __restrict zi) {

                        register __m512 cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcosf(tht);
                        wrkc = _mm512_setzero_ps();
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_ps(cost);
                        csqrt_zmm16r4(divr,divi,&wrkc,&csqr,&csqi);
                        *zr = _mm512_mul_ps(invc,csqr);
                        *zi = _mm512_mul_ps(invc,csqi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                        const float * __restrict __ATTR_ALIGN__(64) pmur,
                                        const float * __restrict __ATTR_ALIGN__(64) pmui,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                        float * __restrict __ATTR_ALIGN__(64) zr,
                                        float * __restrict __ATTR_ALIGN__(64)  zi) {

                        register __m512 tht  = _mm512_load_ps(&ptht[0]);
                        register __m512 mur  = _mm512_load_ps(&pmur[0]);
                        register __m512 mui  = _mm512_load_ps(&pmui[0]);
                        register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512 cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcosf(tht);
                        wrkc = _mm512_setzero_ps();
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_ps(cost);
                        csqrt_zmm16r4(divr,divi,&wrkc,&csqr,&csqi);
                        _mm512_store_ps(&zr[0] ,_mm512_mul_ps(invc,csqr));
                        _mm512_store_ps(&zi[0] ,_mm512_mul_ps(invc,csqi));
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm16r4_u(const float * __restrict  ptht,
                                        const float * __restrict  pmur,
                                        const float * __restrict  pmui,
                                        const float * __restrict  pepsr,
                                        const float * __restrict  pepsi,
                                        float * __restrict  zr,
                                        float * __restrict   zi) {

                        register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                        register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                        register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcosf(tht);
                        wrkc = _mm512_setzero_ps();
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_ps(cost);
                        csqrt_zmm16r4(divr,divi,&wrkc,&csqr,&csqi);
                        _mm512_storeu_ps(&zr[0] ,_mm512_mul_ps(invc,csqr));
                        _mm512_storeu_ps(&zi[0] ,_mm512_mul_ps(invc,csqi));
                }






      } // radiolocation


} // gms


























#endif /*__GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__*/
