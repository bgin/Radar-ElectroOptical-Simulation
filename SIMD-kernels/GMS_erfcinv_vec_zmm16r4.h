

#ifndef __GMS_ERFCINV_VEC_ZMM16R4_H__
#define __GMS_ERFCINV_VEC_ZMM16R4_H__ 081220221614

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

namespace file_version 
{

    const unsigned int GMS_ERFCINV_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_ERFCINV_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_ERFCINV_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_ERFCINV_VEC_ZMM16R4_FULLVER =
      1000U*GMS_ERFCINV_VEC_ZMM16R4_MAJOR+
      100U*GMS_ERFCINV_VEC_ZMM16R4_MINOR+
      10U*GMS_ERFCINV_VEC_ZMM16R4_MICRO;
    const char * const GMS_ERFCINV_VEC_ZMM16R4_CREATION_DATE = "09-12-2022 16:14 AM +00200 (FRI 09 DEC 2022 GMT+2)";
    const char * const GMS_ERFCINV_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_ERFCINV_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_ERFCINV_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized vector of inverse complementary error function values.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"

#if !defined(ERFCINV_VEC_ZMM16R4_PEEL_LOOP)
#define ERFCINV_VEC_ZMM16R4_PEEL_LOOP 1
#endif 

#if !defined(ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH)
#define ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH 1
#endif 

#if !defined(ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS)
#define ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS 1
#endif 

namespace  gms 
{


namespace  math 
{


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void erfcinvv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                   float * __restrict __ATTR_ALIGN__(64) y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   const int32_t n); 
                                                   

               
	         

 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void erfcinvv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                                  float * __restrict  y,
                                                  const __m512 a,
                                                  const __m512 b,
                                                  const __m512 c,
                                                  int32_t n);
              

                 

 

                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void erfcinvv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                   float * __restrict __ATTR_ALIGN__(64) y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   const int32_t n); 

 

                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void erfcinvv_zmm16r4_unroll_6x_u(const float * __restrict  x,
                                                   float * __restrict  y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   int32_t n); 
                

                
	         





              

                    
} // math 
 

} // gms

















#endif /*__GMS_ERFC_VEC_ZMM16R4_H__*/
