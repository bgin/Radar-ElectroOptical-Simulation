

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
			__m512d gust_psd_Ug_zmm8r8(const __m512d sigmau,
			                           const __m512d Lu,
						   const __m512d omega) {

			     __m512d PhiUg;
                             const __m512d invpi    = _mm512_set1_pd(0.318309886183790671537767526745);
			     const __m512d _2       = _mm512_set1_pd(2.0);
			     const __m512d t0       = _mm512_mul_pd(Lu,omega);
			     const __m512d Lomega2  = _mm512_mul_pd(t0,t0);
			     const __m512d _1       = _mm512_set1_pd(1.0);
			     const __m512d sigmau2  = _mm512_mul_pd(sigmau,sigmau);
			     const __m512d L_ov_pi  = _mm512_mul_pd(_mm512_mul_pd(_2,Lu),invpi);
			     const __m512d i_L_omega= _mm512_div_pd(_1,_mm512_add_pd(_1,Lomega));
			     PhiUg                  = _mm512_mul_pd(sigmau2,_mm512_mul_pd(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }



		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512 gust_psd_Ug_zmm16r4(const __m512 sigmau,
			                           const __m512 Lu,
						   const __m512 omega) {

			     __m512 PhiUg;
                             const __m512 invpi    = _mm512_set1_ps(0.318309886183790671537767526745f);
			     const __m512 _2       = _mm512_set1_ps(2.0f);
			     const __m512 t0       = _mm512_mul_ps(Lu,omega);
			     const __m512 Lomega2  = _mm512_mul_ps(t0,t0);
			     const __m512 _1       = _mm512_set1_ps(1.0f);
			     const __m512 sigmau2   = _mm512_mul_ps(sigmau,sigmau);
			     const __m512 L_ov_pi  = _mm512_mul_ps(_mm512_mul_ps(_2,Lu),invpi);
			     const __m512 i_L_omega= _mm512_div_ps(_1,_mm512_add_ps(_1,Lomega));
			     PhiUg                  = _mm512_mul_ps(sigmau2,_mm512_mul_ps(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d gust_psd_Vg_zmm8r8(const __m512d sigmav,
			                           const __m512d Lv,
						   const __m512d omega) {

                           __m512d PhiVg;
			   const __m512d invpi    = _mm512_set1_pd(0.318309886183790671537767526745);
			   const __m512d _2       = _mm512_set1_pd(2.0);
			   const __m512d _12      = _mm512_set1_pd(12.0);
			   const __m512d _4       = _mm512_set1_pd(4.0);
			   const __m512d sigmav2  = _mm512_mul_pd(sigmav,sigmav);
			   const __m512d Lv_o_pi  = _mm512_mul_pd(_mm512_mul_pd(_2,Lv),invpi);
			   const __m512d t0       = _mm512_mul_pd(sigmav2,Lv_o_pi);
			   const __m512d Lvom     = _mm512_mul_pd(Lv,omega);
			   const __m512d Lvom2    = _mm512_mul_pd(Lvom,Lvom);
			   const __m512d num      = _mm512_fmadd_pd(_12,Lvom2,_1);
			   const __m512d denom    = _mm512_fmadd_pd(_4,Lvom2,_1);
			   const __m512d denom2   = _mm512_mul_pd(denom,denom);
			   const __m512d ratio    = _mm512_div_pd(num,denom2);
			   PhiVg                  = _mm512_mul_pd(t0,ratio);
			   return (PhiVg);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512 gust_psd_Vg_zmm16r4(const __m512 sigmav,
			                           const __m512 Lv,
						   const __m512 omega) {

                           __m512 PhiVg;
			   const __m512 invpi    = _mm512_set1_ps(0.318309886183790671537767526745f);
			   const __m512 _2       = _mm512_set1_ps(2.0f);
			   const __m512 _12      = _mm512_set1_ps(12.0f);
			   const __m512 _4       = _mm512_set1_ps(4.0f);
			   const __m512 sigmav2  = _mm512_mul_ps(sigmav,sigmav);
			   const __m512 Lv_o_pi  = _mm512_mul_ps(_mm512_mul_ps(_2,Lv),invpi);
			   const __m512 t0       = _mm512_mul_ps(sigmav2,Lv_o_pi);
			   const __m512 Lvom     = _mm512_mul_ps(Lv,omega);
			   const __m512 Lvom2    = _mm512_mul_ps(Lvom,Lvom);
			   const __m512 num      = _mm512_fmadd_ps(_12,Lvom2,_1);
			   const __m512 denom    = _mm512_fmadd_ps(_4,Lvom2,_1);
			   const __m512 denom2   = _mm512_mul_ps(denom,denom);
			   const __m512 ratio    = _mm512_div_ps(num,denom2);
			   PhiVg                 = _mm512_mul_ps(t0,ratio);
			   return (PhiVg);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d gust_psd_Wg_zmm8r8(const __m512d sigmaw,
			                           const __m512d Lw,
						   const __m512d omega) {
						   
                           __m512d PhiWg = gust_psd_Vg_zmm8r8(sigmaw,Lw,omega);
			   return (PhiWg);
			}


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512 gust_psd_Wg_zmm16r4(const __m512 sigmaw,
			                           const __m512 Lw,
						   const __m512 omega) {
						   
                           __m512 PhiWg = gust_psd_Vg_zmm16r4(sigmaw,Lw,omega);
			   return (PhiWg);
			}



			
    }

}
















#endif /*__GMS_DRYDEN_PSD_AVX512_HPP__*/
