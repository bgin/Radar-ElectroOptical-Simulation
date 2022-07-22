

#ifndef __GMS_BDREF_AVX512_HPP__
#define __GMS_BDREF_AVX512_HPP__ 220720221235

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
      Adapted from the DISORT "BDREF.f" file.
*/

namespace file_version {

    const unsigned int GMS_BDREF_AVX512_MAJOR = 1U;
    const unsigned int GMS_BDREF_AVX512_MINOR = 0U;
    const unsigned int GMS_BDREF_AVX512_MICRO = 0U;
    const unsigned int GMS_BDREF_AVX512_FULLVER =
      1000U*GMS_BDREF_AVX512_MAJOR+
      100U*GMS_BDREF_AVX512_MINOR+
      10U*GMS_BDREF_AVX512_MICRO;
    const char * const GMS_BDREF_AVX512_CREATION_DATE = "22-07-2022 12:35 PM +00200 (FRI 22 JUL 2022 GMT+2)";
    const char * const GMS_BDREF_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_BDREF_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_BDREF_AVX512_DESCRIPTION   = "Vectorized (AVX512) BDREF functions."

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
			__m512d
			bdrf_hapke_zmm8r8(const __m512d mup,
			                  const __m512d mu,
					  const __m512d dphi,
					  const __m512d b0,
					  const __m512d hh,
					  const __m512d w,
					  const __m512d pi) {

                           const __m512d _1   = _mm512_set1_pd(1.0);
			   const __m512d _1_2 = _mm512_set1_pd(0.5);
			   const __m512d _4   = _mm512_set1_pd(4.0);
			   __m512d t0,t1,t2,t3,c0,c1,cdphi,brdf;
			   __m512d calpha,alpha,p,b,h0,gamma,h;
			   cdphi = _mm512_cos_pd(dphi);
			   c0    = _mm512_fmadd_pd(_2,mup,_1);
			   t0 = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                       _mm512_mul_pd(mu,mu)));
			   c1    = _mm512_fmadd_pd(_2,mu,_1); 
			   t1 = _mm512_sqrt_pd(_mm512_sub_pd(_1,
			                       _mm512_mul_pd(mup,mup)));
			   calpha = _mm512_mul_pd(_mm512_fmsub_pd(mu,mup,
			                                      _mm512_mul_pd(t0,t1)),cdphi);
			   gamma  = _mm512_sqrt_pd(_mm512_sub_pd(_1,w);
			   alpha  = _mm512_acos_pd(calpha);
			   p      = _mm512_fmadd_pd(_1_2,calpha,_1);
			   t2     = _mm512_add_pd(hh,_mm512_tan_pd(_mm512_mul_pd(alpha,_1_2)));
			   t3     = _mm512_mul_pd(b0,hh);
			   b      = _mm512_div_pd(t3,t2);
			   h0     = _mm512_div_pd(c0,_mm512_mul_pd(c0,gamma));
			   h      = _mm512_div_pd(c1,_mm512_mul_pd(c1,gamma));
			   t0     = _mm512_div_pd(w,_mm512_mul_pd(_4,pi));
			   t1     = _mm512_add_pd(mu,mup);
			   t2     = _mm512_fmadd_pd(_mm512_add_pd(_1,b),p,
			                                      _mm512_fmsub_pd(h0,h,_1));
			   brdf   = _mm512_mul_pd(_mm512_div_pd(t0,t1),t2);
			   return (brdf);
		      }


		       __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			bdrf_hapke_zmm16r4(const __m512 mup,
			                   const __m512 mu,
					   const __m512 dphi,
					   const __m512 b0,
					   const __m512 hh,
					   const __m512 w,
					   const __m512 pi) {

                           const __m512 _1   = _mm512_set1_ps(1.0f);
			   const __m512 _1_2 = _mm512_set1_ps(0.5f);
			   const __m512 _4   = _mm512_set1_ps(4.0f);
			   __m512 t0,t1,t2,t3,c0,c1,cdphi,brdf;
			   __m512 calpha,alpha,p,b,h0,gamma,h;
			   cdphi = _mm512_cos_ps(dphi);
			   c0    = _mm512_fmadd_ps(_2,mup,_1);
			   t0 = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                       _mm512_mul_ps(mu,mu)));
			   c1    = _mm512_fmadd_ps(_2,mu,_1); 
			   t1 = _mm512_sqrt_ps(_mm512_sub_ps(_1,
			                       _mm512_mul_ps(mup,mup)));
			   calpha = _mm512_mul_ps(_mm512_fmsub_ps(mu,mup,
			                                      _mm512_mul_ps(t0,t1)),cdphi);
			   gamma  = _mm512_sqrt_ps(_mm512_sub_ps(_1,w);
			   alpha  = _mm512_acos_ps(calpha);
			   p      = _mm512_fmadd_ps(_1_2,calpha,_1);
			   t2     = _mm512_add_ps(hh,_mm512_tan_ps(_mm512_mul_ps(alpha,_1_2)));
			   t3     = _mm512_mul_ps(b0,hh);
			   b      = _mm512_div_ps(t3,t2);
			   h0     = _mm512_div_ps(c0,_mm512_mul_ps(c0,gamma));
			   h      = _mm512_div_ps(c1,_mm512_mul_ps(c1,gamma));
			   t0     = _mm512_div_ps(w,_mm512_mul_ps(_4,pi));
			   t1     = _mm512_add_ps(mu,mup);
			   t2     = _mm512_fmadd_ps(_mm512_add_ps(_1,b),p,
			                                      _mm512_fmsub_ps(h0,h,_1));
			   brdf   = _mm512_mul_pd(_mm512_div_pd(t0,t1),t2);
			   return (brdf);
		      }
    }

}


#endif /*__GMS_BDREF_AVX512_HPP__*/
