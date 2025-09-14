

#ifndef __GMS_PHASE_PDF_AVX_H__
#define __GMS_PHASE_PDF_AVX_H__ 121220221332


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
#include "GMS_config.h"



namespace  gms {


           namespace  math {

                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           __ATTR_VECTORCALL__
                   __m256 phase_pdf_ymm8r4(const __m256 a,        
                                           const __m256 phi); 


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           __ATTR_VECTORCALL__
                   __m256d phase_pdf_ymm4r8(const __m256d a,       
                                            const __m256d phi); 


            
         } // math


} // gms













#endif /*__GMS_PHASE_PDF_AVX_H__*/





















