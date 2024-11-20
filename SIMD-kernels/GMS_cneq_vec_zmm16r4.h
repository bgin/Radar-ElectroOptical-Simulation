
#ifndef __GMS_CNEQ_VEC_ZMM16R4_H__
#define __GMS_CNEQ_VEC_ZMM16R4_H__ 221220221843


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

    const unsigned int GMS_CNEQ_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CNEQ_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CNEQ_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CNEQ_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CNEQ_VEC_ZMM16R4_MAJOR+
      100U*GMS_CNEQ_VEC_ZMM16R4_MINOR+
      10U*GMS_CNEQ_VEC_ZMM16R4_MICRO;
    const char * const GMS_CNEQ_VEC_ZMM16R4_CREATION_DATE = "22-12-2022 18:43 AM +00200 (THR 22 DEC 2022 GMT+2)";
    const char * const GMS_CNEQ_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CNEQ_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CNEQ_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex vectors non-equality comparison."

}


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"



namespace  gms {

         namespace math {

                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cneqv_zmm16r4_unroll_16x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  __mmask16 * __restrict eqr,
                                                  __mmask16 * __restrict eqi,
                                                  const int32_t n); 
                                                  
                     

                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cneqv_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                                  __mmask16 * __restrict __ATTR_ALIGN__(64) eqi,
                                                  const int32_t n); 
                                                  

                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cneqv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                                  __mmask16 * __restrict __ATTR_ALIGN__(64) eqi,
                                                  const int32_t n); 

                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cneqv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  __mmask16 * __restrict eqr,
                                                  __mmask16 * __restrict eqi,
                                                  const int32_t n);

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cneqv_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  __mmask16 * __restrict eqr,
                                                  __mmask16 * __restrict eqi,
                                                  const int32_t n); 


                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cneqv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  __mmask16 * __restrict __ATTR_ALIGN__(64) eqr,
                                                  __mmask16 * __restrict __ATTR_ALIGN__(64) eqi,
                                                  const int32_t n); 


      } // math

} // gms












#endif /*__GMS_CNEQ_VEC_ZMM16R4_H__*/
