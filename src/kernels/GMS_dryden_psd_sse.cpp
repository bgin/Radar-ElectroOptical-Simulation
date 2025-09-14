





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


#include <cstdint>
#include "GMS_dryden_psd_sse.h"



			__m128d gms::math::gust_psd_Ug_xmm2r8(const __m128d sigmau,
			                           const __m128d Lu,
						   const __m128d omega) {

			     __m128d PhiUg;
                             const __m128d invpi    = _mm_set1_pd(0.318309886183790671537767526745);
			     const __m128d _2       = _mm_set1_pd(2.0);
			     const __m128d t0       = _mm_mul_pd(Lu,omega);
			     const __m128d Lomega2  = _mm_mul_pd(t0,t0);
			     const __m128d _1       = _mm_set1_pd(1.0);
			     const __m128d sigmau2  = _mm_mul_pd(sigmau,sigmau);
			     const __m128d L_ov_pi  = _mm_mul_pd(_mm_mul_pd(_2,Lu),invpi);
			     const __m128d i_L_omega= _mm_div_pd(_1,_mm_add_pd(_1,Lomega));
			     PhiUg                  = _mm_mul_pd(sigmau2,_mm_mul_pd(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }



		       
			__m128 gms::math::gust_psd_Ug_xmm4r4(const __m128 sigmau,
			                           const __m128 Lu,
						   const __m128 omega) {

			     __m128 PhiUg;
                             const __m128 invpi    = _mm_set1_ps(0.318309886183790671537767526745f);
			     const __m128 _2       = _mm_set1_ps(2.0f);
			     const __m128 t0       = _mm_mul_ps(Lu,omega);
			     const __m128 Lomega2  = _mm_mul_ps(t0,t0);
			     const __m128 _1       = _mm_set1_ps(1.0f);
			     const __m128 sigmau2   = _mm_mul_ps(sigmau,sigmau);
			     const __m128 L_ov_pi  = _mm_mul_ps(_mm_mul_ps(_2,Lu),invpi);
			     const __m128 i_L_omega= _mm_div_ps(_1,_mm_add_ps(_1,Lomega));
			     PhiUg                  = _mm_mul_ps(sigmau2,_mm_mul_ps(L_ov_pi,i_L_omega));
			     return (PhiUg);
		       }


		     
			__m128d gms::math::gust_psd_Vg_xmm2r8(const __m128d sigmav,
			                           const __m128d Lv,
						   const __m128d omega) {

                           __m128d PhiVg;
			   const __m128d invpi    = _mm_set1_pd(0.318309886183790671537767526745);
			   const __m128d _2       = _mm_set1_pd(2.0);
			   const __m128d _12      = _mm_set1_pd(12.0);
			   const __m128d _4       = _mm_set1_pd(4.0);
			   const __m128d sigmav2  = _mm_mul_pd(sigmav,sigmav);
			   const __m128d Lv_o_pi  = _mm_mul_pd(_mm_mul_pd(_2,Lv),invpi);
			   const __m128d t0       = _mm_mul_pd(sigmav2,Lv_o_pi);
			   const __m128d Lvom     = _mm_mul_pd(Lv,omega);
			   const __m128d Lvom2    = _mm_mul_pd(Lvom,Lvom);
			   const __m128d num      = _mm_fmadd_pd(_12,Lvom2,_1);
			   const __m128d denom    = _mm_fmadd_pd(_4,Lvom2,_1);
			   const __m128d denom2   = _mm_mul_pd(denom,denom);
			   const __m128d ratio    = _mm_div_pd(num,denom2);
			   PhiVg                  = _mm_mul_pd(t0,ratio);
			   return (PhiVg);
		       }


		       
			__m128 gms::math::gust_psd_Vg_xmm4r4(const __m128 sigmav,
			                           const __m128 Lv,
						   const __m128 omega) {

                           __m128 PhiVg;
			   const __m128 invpi    = _mm_set1_ps(0.318309886183790671537767526745f);
			   const __m128 _2       = _mm_set1_ps(2.0f);
			   const __m128 _12      = _mm_set1_ps(12.0f);
			   const __m128 _4       = _mm_set1_ps(4.0f);
			   const __m128 sigmav2  = _mm_mul_ps(sigmav,sigmav);
			   const __m128 Lv_o_pi  = _mm_mul_ps(_mm_mul_ps(_2,Lv),invpi);
			   const __m128 t0       = _mm_mul_ps(sigmav2,Lv_o_pi);
			   const __m128 Lvom     = _mm_mul_ps(Lv,omega);
			   const __m128 Lvom2    = _mm_mul_ps(Lvom,Lvom);
			   const __m128 num      = _mm_fmadd_ps(_12,Lvom2,_1);
			   const __m128 denom    = _mm_fmadd_ps(_4,Lvom2,_1);
			   const __m128 denom2   = _mm_mul_ps(denom,denom);
			   const __m128 ratio    = _mm_div_ps(num,denom2);
			   PhiVg                 = _mm_mul_ps(t0,ratio);
			   return (PhiVg);
		       }


		    
			__m128d gms::math::gust_psd_Wg_xmm2r8(const __m128d sigmaw,
			                           const __m128d Lw,
						   const __m128d omega) {
						   
                           __m128d PhiWg = gust_psd_Vg_xmm2r8(sigmaw,Lw,omega);
			   return (PhiWg);
			}


			
			__m128 gms::math::gust_psd_Wg_xmm4r4(const __m128 sigmaw,
			                           const __m128 Lw,
						   const __m128 omega) {
						   
                           __m128 PhiWg = gust_psd_Vg_xmm4r4(sigmaw,Lw,omega);
			   return (PhiWg);
			}




















