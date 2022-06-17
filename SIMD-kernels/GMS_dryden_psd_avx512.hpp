

#ifndef __GMS_DRYDEN_PSD_AVX512_HPP__
#define __GMS_DRYDEN_PSD_AVX512_HPP__ 170620221551



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



/*
    Based on: https://en.wikipedia.org/wiki/Dryden_Wind_Turbulence_Model
*/

namespace file_version {

    const unsigned int GMS_DRYDEN_PSD_AVX512_MAJOR = 1U;
    const unsigned int GMS_DRYDEN_PSD_AVX512_MINOR = 0U;
    const unsigned int GMS_DRYDEN_PSD_AVX512_MICRO = 0U;
    const unsigned int GMS_DRYDEN_PSD_AVX512_FULLVER =
      1000U*GMS_DRYDEN_PSD_AVX512_MAJOR+
      100U*GMS_DRYDEN_PSD_AVX512_MINOR+
      10U*GMS_DRYDEN_PSD_AVX512_MICRO;
    const char * const GMS_DRYDEN_PSD_AVX512_CREATION_DATE = "17-06-2022 15:51 PM +00200 (FRI 17 JUN 2022 GMT+2)";
    const char * const GMS_DRYDEN_PSD_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_DRYDEN_PSD_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_DRYDEN_PSD_AVX512_DESCRIPTION   = "Vectorized (AVX512) Dryden PSD model."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

       namespace math {



                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d gust_psd_Ug_zmm8r8(const __m512d sigma,
			                           const __m512d L,
						   const __m512d omega) {

			     __m512d PhiUg;
                             const __m512d invpi    = _mm512_set1_pd(0.318309886183790671537767526745);
			     const __m512d _2       = _mm512_set1_pd(2.0);
			     const __m512d t0       = _mm512_mul_pd(L,omega);
			     const __m512d Lomega2  = _mm512_mul_pd(t0,t0);
			     const __m512d _1       = _mm512_set1_pd(1.0);
			     const __m512d sigma2   = _mm512_mul_pd(sigma,sigma);
			     const __m512d L_ov_pi  = _mm512_mul_pd(_mm512_mul_pd(_2,L),invpi);
			     const __m512d i_L_omega= _mm512_div_pd(_1,_mm512_add_pd(_1,Lomega));
			     PhiUg                  = _mm512_mul_pd(sigma2,_mm512_mul_pd(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }



		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512 gust_psd_Ug_zmm16r4(const __m512 sigma,
			                           const __m512 L,
						   const __m512 omega) {

			     __m512 PhiUg;
                             const __m512 invpi    = _mm512_set1_ps(0.318309886183790671537767526745f);
			     const __m512 _2       = _mm512_set1_ps(2.0f);
			     const __m512 t0       = _mm512_mul_ps(L,omega);
			     const __m512 Lomega2  = _mm512_mul_ps(t0,t0);
			     const __m512 _1       = _mm512_set1_ps(1.0f);
			     const __m512 sigma2   = _mm512_mul_ps(sigma,sigma);
			     const __m512 L_ov_pi  = _mm512_mul_ps(_mm512_mul_ps(_2,L),invpi);
			     const __m512 i_L_omega= _mm512_div_ps(_1,_mm512_add_ps(_1,Lomega));
			     PhiUg                  = _mm512_mul_ps(sigma2,_mm512_mul_ps(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }


			
    }

}
















#endif /*__GMS_DRYDEN_PSD_AVX512_HPP__*/
