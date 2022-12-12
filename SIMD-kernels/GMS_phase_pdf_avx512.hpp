

#ifndef __GMS_PHASE_PDF_AVX512_HPP__
#define __GMS_PHASE_PDF_AVX512_HPP__ 111220221413


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

    const unsigned int GMS_PHASE_PDF_AVX512_MAJOR = 1U;
    const unsigned int GMS_PHASE_PDF_AVX512_MINOR = 0U;
    const unsigned int GMS_PHASE_PDF_AVX512_MICRO = 0U;
    const unsigned int GMS_PHASE_PDF_AVX512_FULLVER =
      1000U*GMS_PHASE_PDF_AVX512_MAJOR+
      100U*GMS_PHASE_PDF_AVX512_MINOR+
      10U*GMS_PHASE_PDF_AVX512_MICRO;
    const char * const GMS_PHASE_PDF_AVX512_CREATION_DATE = "11-12-2022 14:13 AM +00200 (SUN 11 DEC 2022 GMT+2)";
    const char * const GMS_PHASE_PDF_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_PHASE_PDF_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_PHASE_PDF_AVX512_DESCRIPTION   = "AVX512 optimized Phase PDF."

}

#include <cstdint>
#include <immintrin.h>
#include <cmath>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_sleefsimddp.hpp"
#include "GMS_cephes.h"



namespace  gms {


           namespace  math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   __m512 phase_pdf_zmm16r4(const __m512 a,        // signal-to-noise ratio
                                            const __m512 phi) {    // vector of phi (rad) arguments

                         const __m512 inv2pi = _mm512_set1_ps(0.159154943091895335768883763373f);
                         const __m512 sq2pi  = _mm512_set1_ps(2.506628274631000502415765284811f);
                         const __m512 _1     = _mm512_set1_ps(1.0f);
                         const __m512 _0_5   = _mm512_set1_ps(0.5f);
                         __m512 acphi,errf,expt2,hphi;
                         __m512 aa,haa,invex,t0,t1,t2,c2phi;
                         __m512 f1;
                         aa   = _mm512_mul_ps(a,a);
                         hphi = _mm512_mul_ps(_0_5,phi);
                         haa  = _mm512_mul_ps(_0_5,aa);
                         acphi= _mm512_mul_ps(a,xcosf(phi));
                         t0   = _mm512_mul_ps(sq2pi,acphi);
                         invex= _mm512_div_ps(_1,xexpf(haa));
                         t1   = _mm512_fmadd_ps(sq2pi,acphi,_1);
                         errf = _mm512_erf_ps(acphi);
                         t2   = xcosf(hphi);
                         c2phi= _mm512_mul_ps(t2,t2);
                         expt2= xexpf(_mm512_mul_ps(aa,c2phi));
                         t2   = _mm512_mul_ps(t1,_mm512_mul_ps(errf,expt2));
                         f1   = _mm512_mul_ps(inv2pi,_mm512_mul_ps(invex,t2));
                         return (f1);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   __m512d phase_pdf_zmm8r8(const __m512d a,        // signal-to-noise ratio
                                            const __m512d phi) {    // vector of phi (rad) arguments

                         const __m512d inv2pi = _mm512_set1_pd(0.159154943091895335768883763373);
                         const __m512d sq2pi  = _mm512_set1_pd(2.506628274631000502415765284811);
                         const __m512d _1     = _mm512_set1_pd(1.0f);
                         const __m512d _0_5   = _mm512_set1_pd(0.5f);
                         __m512d acphi,errf,expt2,hphi;
                         __m512d aa,haa,invex,t0,t1,t2,c2phi;
                         __m512d f1;
                         aa   = _mm512_mul_pd(a,a);
                         hphi = _mm512_mul_pd(_0_5,phi);
                         haa  = _mm512_mul_pd(_0_5,aa);
                         acphi= _mm512_mul_pd(a,xcos(phi));
                         t0   = _mm512_mul_pd(sq2pi,acphi);
                         invex= _mm512_div_pd(_1,xexp(haa));
                         t1   = _mm512_fmadd_pd(sq2pi,acphi,_1);
                         errf = _mm512_erf_pd(acphi);
                         t2   = xcos(hphi);
                         c2phi= _mm512_mul_pd(t2,t2);
                         expt2= xexp(_mm512_mul_pd(aa,c2phi));
                         t2   = _mm512_mul_pd(t1,_mm512_mul_pd(errf,expt2));
                         f1   = _mm512_mul_pd(inv2pi,_mm512_mul_pd(invex,t2));
                         return (f1);
                 }


            
         } // math


} // gms













#endif /*__GMS_PHASE_PDF_AVX512_HPP__*/





















