

#ifndef __GMS_PHASE_PDF_SSE_H__
#define __GMS_PHASE_PDF_SSE_H__ 160820240710


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

    const unsigned int GMS_PHASE_PDF_SSE_MAJOR = 1U;
    const unsigned int GMS_PHASE_PDF_SSE_MINOR = 0U;
    const unsigned int GMS_PHASE_PDF_SSE_MICRO = 0U;
    const unsigned int GMS_PHASE_PDF_SSE_FULLVER =
      1000U*GMS_PHASE_PDF_SSE_MAJOR+
      100U*GMS_PHASE_PDF_SSE_MINOR+
      10U*GMS_PHASE_PDF_SSE_MICRO;
    const char * const GMS_PHASE_PDF_SSE_CREATION_DATE = "16-08-2024 07:10 PM +00200 (FRI 16 AUG 2024 07:10PM GMT+2)";
    const char * const GMS_PHASE_PDF_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_PHASE_PDF_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_PHASE_PDF_SSE5_DESCRIPTION   = "SSE optimized Phase PDF."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"



namespace  gms {


           namespace  math {

                   __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m128 phase_pdf_xmm4r4(const __m128 a,       
                                           const __m128 phi); 

                 
	           __ATTR_HOT__
	           __ATTR_VECTORCALL__
                   __m128d phase_pdf_xmm2r8(const __m128d a,        
                                            const __m128d phi); 


            
         } // math


} // gms













#endif /*__GMS_PHASE_PDF_SSE_H__*/





















