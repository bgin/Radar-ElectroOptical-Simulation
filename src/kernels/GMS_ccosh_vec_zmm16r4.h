

#ifndef __GMS_CCOSH_VEC_ZMM16R4_H__
#define __GMS_CCOSH_VEC_ZMM16R4_H__ 271220221214

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

    const unsigned int GMS_CCOSH_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CCOSH_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CCOSH_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CCOSH_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CCOSH_VEC_ZMM16R4_MAJOR+
      100U*GMS_CCOSH_VEC_ZMM16R4_MINOR+
      10U*GMS_CCOSH_VEC_ZMM16R4_MICRO;
    const char * const GMS_CCOSH_VEC_ZMM16R4_CREATION_DATE = "27-12-2022 12:14  +00200 (TUE 27 DEC 2022 GMT+2)";
    const char * const GMS_CCOSH_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CCOSH_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CCOSH_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized complex cosine hyperbolic function."

}

#include <cstdint>
#include "GMS_config.h"




namespace  gms {


      namespace math {

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	            void ccoshv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict csre,
                                                   float * __restrict csim,
                                                   const int32_t n); 


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void ccoshv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) csre,
                                                   float * __restrict __ATTR_ALIGN__(64) csim,
                                                   const int32_t n); 

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void ccoshv_zmm16r4_unroll_8x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict csre,
                                                   float * __restrict csim,
                                                   const int32_t n); 

                  
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void ccoshv_zmm16r4_unroll_8x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) csre,
                                                   float * __restrict __ATTR_ALIGN__(64) csim,
                                                   const int32_t n); 
                  
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void ccoshv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict  __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) csre,
                                                   float * __restrict __ATTR_ALIGN__(64) csim,
                                                   const int32_t n); 

                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void ccoshv_zmm16r4_unroll_8x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict csre,
                                                   float * __restrict csim,
                                                   const int32_t n); 
                 


                 
                  

     } // math

} // gms





#endif /*__GMS_CCOSH_VEC_ZMM16R4_H__*/
