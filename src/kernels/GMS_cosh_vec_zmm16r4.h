

#ifndef __GMS_COSH_VEC_ZMM16R4_H__
#define __GMS_COSH_VEC_ZMM16R4_H__ 090120220833

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

    const unsigned int GMS_COSH_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_COSH_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_COSH_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_COSH_VEC_ZMM16R4_FULLVER =
      1000U*GMS_COSH_VEC_ZMM16R4_MAJOR+
      100U*GMS_COSH_VEC_ZMM16R4_MINOR+
      10U*GMS_COSH_VEC_ZMM16R4_MICRO;
    const char * const GMS_COSH_VEC_ZMM16R4_CREATION_DATE = "09-12-2022 08:33 AM +00200 (FRI 09 DEC 2022 GMT+2)";
    const char * const GMS_COSH_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COSH_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COSH_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized vector of cosine hyperbolic values."

}

#include <cstdint>
#include "GMS_config.h"


namespace  gms {


         namespace  math {

/*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                   float * __restrict __ATTR_ALIGN__(64) y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   const int32_t n); 



                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_mask_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x, 
                                                       const float * __restrict __ATTR_ALIGN__(64) z,
                                                       float * __restrict __ATTR_ALIGN__(64) y,
                                                       const __mmask16 * __restrict __ATTR_ALIGN__(64) m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive); 
                                                       

              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_mask_zmm16r4_unroll_10x_u(const float * __restrict  x, 
                                                       const float * __restrict  z,
                                                       float * __restrict  y,
                                                       const __mmask16 *  m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive); 


              /*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
             
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                                  float * __restrict  y,
                                                  const __m512 a,
                                                  const __m512 b,
                                                  const __m512 c,
                                                  const int32_t n); 
               
              

                 

 /*             
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/

               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                   float * __restrict __ATTR_ALIGN__(64) y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   const int32_t n);
               /*             
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_zmm16r4_unroll_6x_u(const float * __restrict  x,
                                                   float * __restrict  y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   const int32_t n); 

                

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_mask_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x, 
                                                       const float * __restrict __ATTR_ALIGN__(64) z,
                                                       float * __restrict __ATTR_ALIGN__(64) y,
                                                       const __mmask16 * __restrict __ATTR_ALIGN__(64) m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive);

                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void coshv_mask_zmm16r4_unroll_6x_u(const float * __restrict  x, 
                                                       const float * __restrict  z,
                                                       float * __restrict  y,
                                                       const __mmask16 *  m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive); 





              

                    

        } // math 
 

} // gms

















#endif /*__GMS_COSH_VEC_ZMM16R4_H__*/
