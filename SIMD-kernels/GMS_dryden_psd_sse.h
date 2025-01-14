

#ifndef __GMS_DRYDEN_PSD_SSE_H__
#define __GMS_DRYDEN_PSD_SSE_H__ 011120230855



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

    const unsigned int GMS_DRYDEN_PSD_SSE_MAJOR = 1U;
    const unsigned int GMS_DRYDEN_PSD_SSE_MINOR = 0U;
    const unsigned int GMS_DRYDEN_PSD_SSE_MICRO = 0U;
    const unsigned int GMS_DRYDEN_PSD_SSE_FULLVER =
      1000U*GMS_DRYDEN_PSD_SSE_MAJOR+
      100U*GMS_DRYDEN_PSD_SSE_MINOR+
      10U*GMS_DRYDEN_PSD_SSE_MICRO;
    const char * const GMS_DRYDEN_PSD_SSE_CREATION_DATE = "01-11-2023 08:55 PM +00200 (WED 01 NOV 2023 GMT+2)";
    const char * const GMS_DRYDEN_PSD_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_DRYDEN_PSD_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_DRYDEN_PSD_SSE_DESCRIPTION   = "Vectorized (SSE) Dryden PSD model."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

       namespace math {



                    
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m128d gust_psd_Ug_xmm2r8(const __m128d sigmau,
			                           const __m128d Lu,
						   const __m128d omega); 



		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m128 gust_psd_Ug_xmm4r4(const __m128 sigmau,
			                           const __m128 Lu,
						   const __m128 omega); 


		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m128d gust_psd_Vg_xmm2r8(const __m128d sigmav,
			                           const __m128d Lv,
						   const __m128d omega); 


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m128 gust_psd_Vg_xmm4r4(const __m128 sigmav,
			                           const __m128 Lv,
						   const __m128 omega); 


		     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m128d gust_psd_Wg_xmm2r8(const __m128d sigmaw,
			                           const __m128d Lw,
						   const __m128d omega); 

			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m128 gust_psd_Wg_xmm4r4(const __m128 sigmaw,
			                           const __m128 Lw,
						   const __m128 omega); 


			
    }

}
















#endif /*__GMS_DRYDEN_PSD_SSE_H__*/
