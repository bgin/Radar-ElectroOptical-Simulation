

#ifndef __GMS_CDIV_VEC_ZMM16R4_H__
#define __GMS_CDIV_VEC_ZMM16R4_H__


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

    const unsigned int GMS_CDIV_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_CDIV_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_CDIV_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_CDIV_VEC_ZMM16R4_FULLVER =
      1000U*GMS_CDIV_VEC_ZMM16R4_MAJOR+
      100U*GMS_CDIV_VEC_ZMM16R4_MINOR+
      10U*GMS_CDIV_VEC_ZMM16R4_MICRO;
    const char * const GMS_CDIV_VEC_ZMM16R4_CREATION_DATE = "15-12-2022 13:05 AM +00200 (THR 15 DEC 2022 GMT+2)";
    const char * const GMS_CDIV_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_CDIV_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_CDIV_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized decomposed complex vector division operations."

}

#include <cstdint>
#include "GMS_config.h"



namespace  gms {

          namespace math {

                  
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdivv_zmm16r4_unroll_10x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n); 
              ///////////////////////////////////////////////////////////////////////////////////
              
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdivv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict       __ATTR_ALIGN__(64) zre,
                                                  float * __restrict       __ATTR_ALIGN__(64) zim,
                                                  const int32_t n);

                 /////////////////////////////////////////////////////////////////////////////
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdivv_zmm16r4_unroll_8x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n); 

              ///////////////////////////////////////////////////////////////////////////////////

                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdivv_zmm16r4_unroll_8x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict       __ATTR_ALIGN__(64) zre,
                                                  float * __restrict       __ATTR_ALIGN__(64) zim,
                                                  const int32_t n); 
                         
 
                ////////////////////////////////////////////////////////////////////////// 
                     
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdivv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict       __ATTR_ALIGN__(64) zre,
                                                  float * __restrict       __ATTR_ALIGN__(64) zim,
                                                  const int32_t n); 
                /////////////////////////////////////////////////////////////////////

                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdivv_zmm16r4_unroll_6x_u(const float * __restrict  xre,
                                                  const float * __restrict  xim,
                                                  const float * __restrict  yre,
                                                  const float * __restrict  yim,
                                                  float * __restrict        zre,
                                                  float * __restrict        zim,
                                                  const int32_t n); 
                     
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdiv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n); 

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	          void cdiv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n); 

                  
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdiv_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                  const float * __restrict xim, 
                                                  const __m512             vs,
                                                  float  * __restrict zre,
                                                  float  * __restrict zim,
                                                  const int32_t n);

                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void cdiv_zmm16r4_unroll_6x_a( const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim, 
                                                  const __m512             vs,
                                                  float  * __restrict __ATTR_ALIGN__(64) zre,
                                                  float  * __restrict __ATTR_ALIGN__(64) zim,
                                                  const int32_t n); 



       }

}












#endif /*__GMS_CDIV_VEC_ZMM16R4_H__*/
