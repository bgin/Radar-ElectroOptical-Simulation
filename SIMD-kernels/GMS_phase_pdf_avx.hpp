

#ifndef __GMS_PHASE_PDF_AVX_HPP__
#define __GMS_PHASE_PDF_AVX_HPP__ 121220221332


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

    const unsigned int GMS_PHASE_PDF_AVX_MAJOR = 1U;
    const unsigned int GMS_PHASE_PDF_AVX_MINOR = 0U;
    const unsigned int GMS_PHASE_PDF_AVX_MICRO = 0U;
    const unsigned int GMS_PHASE_PDF_AVX_FULLVER =
      1000U*GMS_PHASE_PDF_AVX_MAJOR+
      100U*GMS_PHASE_PDF_AVX_MINOR+
      10U*GMS_PHASE_PDF_AVX_MICRO;
    const char * const GMS_PHASE_PDF_AVX_CREATION_DATE = "12-12-2022 13:32 AM +00200 (MON 12 DEC 2022 GMT+2)";
    const char * const GMS_PHASE_PDF_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_PHASE_PDF_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_PHASE_PDF_AVX_DESCRIPTION   = "AVX optimized Phase PDF."

}

#include <cstdint>
#include <immintrin.h>




namespace  gms {


           namespace  math {

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   __m256 phase_pdf_ymm8r4(const __m256 a,        // signal-to-noise ratio
                                           const __m256 phi) {    // vector of phi (rad) arguments

                         const __m256 inv2pi = _mm256_set1_ps(0.159154943091895335768883763373f);
                         const __m256 sq2pi  = _mm256_set1_ps(2.506628274631000502415765284811f);
                         const __m256 _1     = _mm256_set1_ps(1.0f);
                         const __m256 _0_5   = _mm256_set1_ps(0.5f);
                         __m256 acphi,errf,expt2,hphi;
                         __m256 aa,haa,invex,t0,t1,t2,c2phi;
                         __m256 f1;
                         aa   = _mm256_mul_ps(a,a);
                         hphi = _mm256_mul_ps(_0_5,phi);
                         haa  = _mm256_mul_ps(_0_5,aa);
                         acphi= _mm256_mul_ps(a,_mm256_cos_ps(phi));
                         t0   = _mm256_mul_ps(sq2pi,acphi);
                         invex= _mm256_div_ps(_1,_mm256_exp_ps(haa));
                         t1   = _mm256_fmadd_ps(sq2pi,acphi,_1);
                         errf = _mm256_erf_ps(acphi);
                         t2   = _mm256_cos_ps(hphi);
                         c2phi= _mm256_mul_ps(t2,t2);
                         expt2= _mm256_exp_ps(_mm256_mul_ps(aa,c2phi));
                         t2   = _mm256_mul_ps(t1,_mm256_mul_ps(errf,expt2));
                         f1   = _mm256_mul_ps(inv2pi,_mm256_mul_ps(invex,t2));
                         return (f1);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   __m256d phase_pdf_ymm4r8(const __m256d a,        // signal-to-noise ratio
                                            const __m256d phi) {    // vector of phi (rad) arguments

                         const __m256d inv2pi = _mm256_set1_pd(0.159154943091895335768883763373);
                         const __m256d sq2pi  = _mm256_set1_pd(2.506628274631000502415765284811);
                         const __m256d _1     = _mm256_set1_pd(1.0f);
                         const __m256d _0_5   = _mm256_set1_pd(0.5f);
                         __m256d acphi,errf,expt2,hphi;
                         __m256d aa,haa,invex,t0,t1,t2,c2phi;
                         __m256d f1;
                         aa   = _mm256_mul_pd(a,a);
                         hphi = _mm256_mul_pd(_0_5,phi);
                         haa  = _mm256_mul_pd(_0_5,aa);
                         acphi= _mm256_mul_pd(a,_mm256_cos_pd(phi));
                         t0   = _mm256_mul_pd(sq2pi,acphi);
                         invex= _mm256_div_pd(_1,_mm256_exp_pd(haa));
                         t1   = _mm256_fmadd_pd(sq2pi,acphi,_1);
                         errf = _mm256_erf_pd(acphi);
                         t2   = _mm256_cos_pd(hphi);
                         c2phi= _mm256_mul_pd(t2,t2);
                         expt2= _mm256_exp_pd(_mm256_mul_pd(aa,c2phi));
                         t2   = _mm256_mul_pd(t1,_mm256_mul_pd(errf,expt2));
                         f1   = _mm256_mul_pd(inv2pi,_mm256_mul_pd(invex,t2));
                         return (f1);
                 }


            
         } // math


} // gms













#endif /*__GMS_PHASE_PDF_AVX_HPP__*/





















