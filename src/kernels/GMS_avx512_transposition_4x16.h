

#ifndef __GMS_AVX512_TRANSPOSITION_4X16_HPP__
#define __GMS_AVX512_TRANSPOSITION_4X16_HPP__ 240420230952

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

    const unsigned int GMS_AVX512_TRANSPOSITION_4X16_MAJOR = 1U;
    const unsigned int GMS_AVX512_TRANSPOSITION_4X16_MINOR = 0U;
    const unsigned int GMS_AVX512_TRANSPOSITION_4X16_MICRO = 0U;
    const unsigned int GMS_AVX512_TRANSPOSITION_4X16_FULLVER =
      1000U*GMS_AVX512_TRANSPOSITION_4X16_MAJOR+
      100U*GMS_AVX512_TRANSPOSITION_4X16_MINOR+
      10U*GMS_AVX512_TRANSPOSITION_4X16_MICRO;
    const char * const GMS_AVX512_TRANSPOSITION_4X16_CREATION_DATE = "24-04-2023 09:52 AM +00200 (MON 24 APR 2023 GMT+2)";
    const char * const GMS_AVX512_TRANSPOSITION_4X16_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_AVX512_TRANSPOSITION_4X16_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_AVX512_TRANSPOSITION_4X16_DESCRIPTION   = "AVX512 optimized matrix 24x16 transposition kernels."

}


#include <immintrin.h>
#include "GMS_config.h"


namespace gms {


         namespace math {

 
                     /*
                            In-place.
                            Data pointer must point to an array (__m512) of 4 elements.
                       */

                                        
		     
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16(__m512 * __restrict inout);

                     /*
                            Out-of-place.
                            Pointers must point to an arrays (__m512) of 4 elements.
                       */

                     
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16(const __m512 * __restrict in,
                                                  __m512 * __restrict out);
                                                  

		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16_u(float * __restrict inout);


                   
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16_a(float * __restrict __ATTR_ALIGN__(64) inout);


                     
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16_u(float * __restrict in,
                                                    float * __restrict out);


                      
                     
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16_a(float * __restrict __ATTR_ALIGN__(64) in,
                                                    float * __restrict __ATTR_ALIGN__(64) out);


                    
		       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      __ATTR_VECTORCALL__
		      void transpose_zmm16r4_4x16(const __m512 in0,
                                                  const __m512 in1,
                                                  const __m512 in2,
                                                  const __m512 in3,
                                                  __m512 * __restrict out0,
                                                  __m512 * __restrict out1,
                                                  __m512 * __restrict out2,
                                                  __m512 * __restrict out3);

                    

#include <cstdint>


                     
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void transpose_zmm16r4_4x16_u8x_u(float * __restrict in,
                                                        float * __restrict out,
                                                        const int32_t n,
                                                        const int32_t m);
                                                           

		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void transpose_zmm16r4_4x16_u4x_u(float * __restrict in,
                                                        float * __restrict out,
                                                        const int32_t n,
                                                        const int32_t m); 
                                                        

                    
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void transpose_zmm16r4_4x16_u8x_a(float * __restrict __ATTR_ALIGN__(64) in,
                                                        float * __restrict __ATTR_ALIGN__(64) out,
                                                        const int32_t n,
                                                        const int32_t m); 

                           

                    
		      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void transpose_zmm16r4_4x16_u4x_a(float * __restrict __ATTR_ALIGN__(64) in,
                                                        float * __restrict __ATTR_ALIGN__(64) out,
                                                        const int32_t n,
                                                        const int32_t m); 
                     
                                                      

       } // math




} // gms























#endif /*__GMS_AVX512_TRANSPOSITION_4X16_H__*/
