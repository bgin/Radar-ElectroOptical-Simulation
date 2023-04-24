

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

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_zmm16r4_4x16(__m512 * __restrict inout) {

                           _mm_prefetch((const char*)&inout[0],_MM_HINT_T0);
                           register __m512 x0 = _mm512_unpacklo_ps(inout[0],inout[1]);
                           register __m512 x1 = _mm512_unpackhi_ps(inout[0],inout[1]);
                           _mm_prefetch((const char*)&inout[2],_MM_HINT_T0);
                           register __m512 x2 = _mm512_unpacklo_ps(inout[2],inout[3]);
                           register __m512 x3 = _mm512_unpackhi_ps(inout[2],inout[3]);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           inout[0]        = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0));
                           inout[1]        = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0));
                           inout[2]        = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1));
                           inout[3]        = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1));
                    }

                     /*
                            Out-of-place.
                            Pointers must point to an arrays (__m512) of 4 elements.
                       */

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_zmm16r4_4x16(const __m512 * __restrict in,
                                                  __m512 * __restrict out) {

                           _mm_prefetch((const char*)&in[0],_MM_HINT_T0);
                           register __m512 x0 = _mm512_unpacklo_ps(in[0],in[1]);
                           register __m512 x1 = _mm512_unpackhi_ps(in[0],in[1]);
                           _mm_prefetch((const char*)&in[2],_MM_HINT_T0);
                           register __m512 x2 = _mm512_unpacklo_ps(in[2],in[3]);
                           register __m512 x3 = _mm512_unpackhi_ps(in[2],in[3]);
                           register __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           register __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           register __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           out[0]         = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0));
                           out[1]         = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0));
                           out[2]         = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1));
                           out[3]         = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1));
                    }


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void transpose_zmm16r4_4x16_u(float * __restrict inout) {

                           _mm_prefetch((const char*)&inout[0],_MM_HINT_T0);
                            __
                            __m512 x0 = _mm512_unpacklo_ps(inout[0],inout[1]);
                            __m512 x1 = _mm512_unpackhi_ps(inout[0],inout[1]);
                           _mm_prefetch((const char*)&inout[2],_MM_HINT_T0);
                           __m512 x2 = _mm512_unpacklo_ps(inout[2],inout[3]);
                           __m512 x3 = _mm512_unpackhi_ps(inout[2],inout[3]);
                           __m512 y0       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(1,0,1,0));
                           __m512 y1       = _mm512_shuffle_ps(x0,x2,_MM_SHUFFLE(3,2,3,2));
                           __m512 y2       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(1,0,1,0));
                           __m512 y3       = _mm512_shuffle_ps(x1,x3,_MM_SHUFFLE(3,2,3,2));
                                  x0       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(2,0,2,0));
                                  x1       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(2,0,2,0));
                                  x2       = _mm512_shuffle_f32x4(y0,y1,_MM_SHUFFLE(3,1,3,1));
                                  x3       = _mm512_shuffle_f32x4(y2,y3,_MM_SHUFFLE(3,1,3,1));
                           inout[0]        = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(2,0,2,0));
                           inout[1]        = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(2,0,2,0));
                           inout[2]        = _mm512_shuffle_f32x4(x0,x1,_MM_SHUFFLE(3,1,3,1));
                           inout[3]        = _mm512_shuffle_f32x4(x2,x3,_MM_SHUFFLE(3,1,3,1));
                    }


                     
                                                      

       } // math




} // gms























#endif /*__GMS_AVX512_TRANSPOSITION_4X16_HPP__*/
