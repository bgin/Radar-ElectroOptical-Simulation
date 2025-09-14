

#ifndef __GMS_BDREF_AVX_H__
#define __GMS_BDREF_AVX_H__ 

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

    const unsigned int GMS_BDREF_AVX_MAJOR = 1U;
    const unsigned int GMS_BDREF_AVX_MINOR = 0U;
    const unsigned int GMS_BDREF_AVX_MICRO = 0U;
    const unsigned int GMS_BDREF_AVX_FULLVER =
      1000U*GMS_BDREF_AVX_MAJOR+
      100U*GMS_BDREF_AVX_MINOR+
      10U*GMS_BDREF_AVX_MICRO;
    const char * const GMS_BDREF_AVX_CREATION_DATE = "30-10-2023 15:18 +00200 (MON 30 OCT 2023 GMT+2)";
    const char * const GMS_BDREF_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_BDREF_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_BDREF_AVX_DESCRIPTION   = "Vectorized (AVX/AVX2) BDREF functions."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

      namespace math {

                       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	              	__m256d
			bdrf_hapke_ymm4r8(const __m256d mup,
			                  const __m256d mu,
					  const __m256d dphi,
					  const __m256d b0,
					  const __m256d hh,
					  const __m256d w,
					  const __m256d pi); 
					  
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m256
			bdrf_hapke_ymm8r4(const __m256 mup,
			                   const __m256 mu,
					   const __m256 dphi,
					   const __m256 b0,
					   const __m256 hh,
					   const __m256 w,
					   const __m256 pi); 
					   
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                 __m256d
			bdrf_rpv_ymm4r8(const __m256d mu_i,
			                const __m256d mu_r,
					const __m256d dphi,
					const __m256d rh00,
					const __m256d kappa,
					const __m256d g_hg,
					const __m256d h0); 
					
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m256
			bdrf_rpv_ymm8r4(const __m256 mu_i,
			                const __m256 mu_r,
					const __m256 dphi,
					const __m256 rh00,
					const __m256 kappa,
					const __m256 g_hg,
					const __m256 h0); 


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	             	__m256d
			bdrf_rossli_ymm4r8(const __m256d mu_i,
			                   const __m256d mu_r,
					   const __m256d dphi,
					   const __m256d k_iso,
					   const __m256d k_vol,
					   const __m256d k_geo,
					   const __m256d alpha0); 
					   
					   
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                __m256
			bdrf_rossli_ymm8r4(const __m256 mu_i,
			                   const __m256 mu_r,
					   const __m256 dphi,
					   const __m256 k_iso,
					   const __m256 k_vol,
					   const __m256 k_geo,
					   const __m256 alpha0); 
					   
					   
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	              	__m256d
			brdf_ocean_ymm4r8(const bool do_shadow,
			                  const __m256d refrac_idx,
					  const __m256d ws,
					  const __m256d mu_i,
					  const __m256d mu_r,
					  const __m256d dphi); 
					  
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m256d
			shadow_eta_ymm4r8(const __m256d cos_theta,
			                  const __m256d sigma_sq,
					  const __m256d pi); 
					  

		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m256
			brdf_ocean_ymm8r4(const bool do_shadow,
			                  const __m256 refrac_idx,
					  const __m256 ws,
					  const __m256 mu_i,
					  const __m256 mu_r,
					  const __m256 dphi); 
					  


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               	__m256
			shadow_eta_ymm8r4(const __m256 cos_theta,
			                  const __m256 sigma_sq,
					  const __m256 pi);
					 




		      
    }

}


#endif /*__GMS_BDREF_AVX_HPP__*/
