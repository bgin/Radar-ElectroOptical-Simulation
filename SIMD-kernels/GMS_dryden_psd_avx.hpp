

#ifndef __GMS_DRYDEN_PSD_AVX_HPP__
#define __GMS_DRYDEN_PSD_AVX_HPP__ 011120230849



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

    const unsigned int GMS_DRYDEN_PSD_AVX_MAJOR = 1U;
    const unsigned int GMS_DRYDEN_PSD_AVX_MINOR = 0U;
    const unsigned int GMS_DRYDEN_PSD_AVX_MICRO = 0U;
    const unsigned int GMS_DRYDEN_PSD_AVX_FULLVER =
      1000U*GMS_DRYDEN_PSD_AVX_MAJOR+
      100U*GMS_DRYDEN_PSD_AVX_MINOR+
      10U*GMS_DRYDEN_PSD_AVX_MICRO;
    const char * const GMS_DRYDEN_PSD_AVX_CREATION_DATE = "01-11-2023 08:49 PM +00200 (WED 01 NOV 2023 GMT+2)";
    const char * const GMS_DRYDEN_PSD_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_DRYDEN_PSD_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_DRYDEN_PSD_AVX_DESCRIPTION   = "Vectorized (AVX) Dryden PSD model."

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
			__m256d gust_psd_Ug_ymm4r8(const __m256d sigmau,
			                           const __m256d Lu,
						   const __m256d omega) {

			     __m256d PhiUg;
                             const __m256d invpi    = _mm256_set1_pd(0.318309886183790671537767526745);
			     const __m256d _2       = _mm256_set1_pd(2.0);
			     const __m256d t0       = _mm256_mul_pd(Lu,omega);
			     const __m256d Lomega2  = _mm256_mul_pd(t0,t0);
			     const __m256d _1       = _mm256_set1_pd(1.0);
			     const __m256d sigmau2  = _mm256_mul_pd(sigmau,sigmau);
			     const __m256d L_ov_pi  = _mm256_mul_pd(_mm256_mul_pd(_2,Lu),invpi);
			     const __m256d i_L_omega= _mm256_div_pd(_1,_mm256_add_pd(_1,Lomega));
			     PhiUg                  = _mm256_mul_pd(sigmau2,_mm256_mul_pd(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }



		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256 gust_psd_Ug_ymm8r4(const __m256 sigmau,
			                           const __m256 Lu,
						   const __m256 omega) {

			     __m256 PhiUg;
                             const __m256 invpi    = _mm256_set1_ps(0.318309886183790671537767526745f);
			     const __m256 _2       = _mm256_set1_ps(2.0f);
			     const __m256 t0       = _mm256_mul_ps(Lu,omega);
			     const __m256 Lomega2  = _mm256_mul_ps(t0,t0);
			     const __m256 _1       = _mm256_set1_ps(1.0f);
			     const __m256 sigmau2   = _mm256_mul_ps(sigmau,sigmau);
			     const __m256 L_ov_pi  = _mm256_mul_ps(_mm256_mul_ps(_2,Lu),invpi);
			     const __m256 i_L_omega= _mm256_div_ps(_1,_mm256_add_ps(_1,Lomega));
			     PhiUg                  = _mm256_mul_ps(sigmau2,_mm256_mul_ps(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256d gust_psd_Vg_ymm4r8(const __m256d sigmav,
			                           const __m256d Lv,
						   const __m256d omega) {

                           __m256d PhiVg;
			   const __m256d invpi    = _mm256_set1_pd(0.318309886183790671537767526745);
			   const __m256d _2       = _mm256_set1_pd(2.0);
			   const __m256d _12      = _mm256_set1_pd(12.0);
			   const __m256d _4       = _mm256_set1_pd(4.0);
			   const __m256d sigmav2  = _mm256_mul_pd(sigmav,sigmav);
			   const __m256d Lv_o_pi  = _mm256_mul_pd(_mm256_mul_pd(_2,Lv),invpi);
			   const __m256d t0       = _mm256_mul_pd(sigmav2,Lv_o_pi);
			   const __m256d Lvom     = _mm256_mul_pd(Lv,omega);
			   const __m256d Lvom2    = _mm256_mul_pd(Lvom,Lvom);
			   const __m256d num      = _mm256_fmadd_pd(_12,Lvom2,_1);
			   const __m256d denom    = _mm256_fmadd_pd(_4,Lvom2,_1);
			   const __m256d denom2   = _mm256_mul_pd(denom,denom);
			   const __m256d ratio    = _mm256_div_pd(num,denom2);
			   PhiVg                  = _mm256_mul_pd(t0,ratio);
			   return (PhiVg);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256 gust_psd_Vg_ymm8r4(const __m256 sigmav,
			                           const __m256 Lv,
						   const __m256 omega) {

                           __m256 PhiVg;
			   const __m256 invpi    = _mm256_set1_ps(0.318309886183790671537767526745f);
			   const __m256 _2       = _mm256_set1_ps(2.0f);
			   const __m256 _12      = _mm256_set1_ps(12.0f);
			   const __m256 _4       = _mm256_set1_ps(4.0f);
			   const __m256 sigmav2  = _mm256_mul_ps(sigmav,sigmav);
			   const __m256 Lv_o_pi  = _mm256_mul_ps(_mm256_mul_ps(_2,Lv),invpi);
			   const __m256 t0       = _mm256_mul_ps(sigmav2,Lv_o_pi);
			   const __m256 Lvom     = _mm256_mul_ps(Lv,omega);
			   const __m256 Lvom2    = _mm256_mul_ps(Lvom,Lvom);
			   const __m256 num      = _mm256_fmadd_ps(_12,Lvom2,_1);
			   const __m256 denom    = _mm256_fmadd_ps(_4,Lvom2,_1);
			   const __m256 denom2   = _mm256_mul_ps(denom,denom);
			   const __m256 ratio    = _mm256_div_ps(num,denom2);
			   PhiVg                 = _mm256_mul_ps(t0,ratio);
			   return (PhiVg);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256d gust_psd_Wg_ymm4r8(const __m256d sigmaw,
			                           const __m256d Lw,
						   const __m256d omega) {
						   
                           __m256d PhiWg = gust_psd_Vg_ymm4r8(sigmaw,Lw,omega);
			   return (PhiWg);
			}


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256 gust_psd_Wg_ymm8r4(const __m256 sigmaw,
			                           const __m256 Lw,
						   const __m256 omega) {
						   
                           __m256 PhiWg = gust_psd_Vg_ymm8r4(sigmaw,Lw,omega);
			   return (PhiWg);
			}



			
    }

}
















#endif /*__GMS_DRYDEN_PSD_AVX_HPP__*/
