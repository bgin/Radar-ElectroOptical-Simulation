

#ifndef __GMS_BDREF_SSE_H__
#define __GMS_BDREF_SSE_H__ 291020231518

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

    const unsigned int GMS_BDREF_SSE_MAJOR = 1U;
    const unsigned int GMS_BDREF_SSE_MINOR = 0U;
    const unsigned int GMS_BDREF_SSE_MICRO = 0U;
    const unsigned int GMS_BDREF_SSE_FULLVER =
      1000U*GMS_BDREF_SSE_MAJOR+
      100U*GMS_BDREF_SSE_MINOR+
      10U*GMS_BDREF_SSE_MICRO;
    const char * const GMS_BDREF_SSE_CREATION_DATE = "29-10-2023 15:18 +00200 (SUN 29 OCT 2023 GMT+2)";
    const char * const GMS_BDREF_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_BDREF_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_BDREF_SSE_DESCRIPTION   = "Vectorized (SSE) BDREF functions."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

      namespace math {

                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                __m128d
			bdrf_hapke_xmm2r8(const __m128d mup,
			                  const __m128d mu,
					  const __m128d dphi,
					  const __m128d b0,
					  const __m128d hh,
					  const __m128d w,
					  const __m128d pi); 


		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	              	__m128
			bdrf_hapke_xmm4r4(const __m128 mup,
			                   const __m128 mu,
					   const __m128 dphi,
					   const __m128 b0,
					   const __m128 hh,
					   const __m128 w,
					   const __m128 pi); 


                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m128d
			bdrf_rpv_xmm2r8(const __m128d mu_i,
			                const __m128d mu_r,
					const __m128d dphi,
					const __m128d rh00,
					const __m128d kappa,
					const __m128d g_hg,
					const __m128d h0); 


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m128
			bdrf_rpv_xmm4r4(const __m128 mu_i,
			                const __m128 mu_r,
					const __m128 dphi,
					const __m128 rh00,
					const __m128 kappa,
					const __m128 g_hg,
					const __m128 h0);
					


		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m128d
			bdrf_rossli_xmm2r8(const __m128d mu_i,
			                   const __m128d mu_r,
					   const __m128d dphi,
					   const __m128d k_iso,
					   const __m128d k_vol,
					   const __m128d k_geo,
					   const __m128d alpha0); 
					   

		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	              	__m128
			bdrf_rossli_xmm4r4(const __m128 mu_i,
			                   const __m128 mu_r,
					   const __m128 dphi,
					   const __m128 k_iso,
					   const __m128 k_vol,
					   const __m128 k_geo,
					   const __m128 alpha0);
					   


		     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m128d
			brdf_ocean_xmm2r8(const bool do_shadow,
			                  const __m128d refrac_idx,
					  const __m128d ws,
					  const __m128d mu_i,
					  const __m128d mu_r,
					  const __m128d dphi);
					  

		     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m128d
			shadow_eta_xmm2r8(const __m128d cos_theta,
			                  const __m128d sigma_sq,
					  const __m128d pi); 
					 


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                __m128
			brdf_ocean_xmm4r4(const bool do_shadow,
			                  const __m128 refrac_idx,
					  const __m128 ws,
					  const __m128 mu_i,
					  const __m128 mu_r,
					  const __m128 dphi); 
					  

		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                __m128
			shadow_eta_xmm4r4(const __m128 cos_theta,
			                  const __m128 sigma_sq,
					  const __m128 pi); 




		      
    }

}


#endif /*__GMS_BDREF_SSE_H__*/
