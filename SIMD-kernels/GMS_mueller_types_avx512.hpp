

#ifndef __GMS_MUELLER_TYPES_AVX512_HPP__
#define __GMS_MUELLER_TYPES_AVX512_HPP__ 211220211503
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


namespace file_info {

 const unsigned int gGMS_MUELLER_TYPES_AVX512_MAJOR = 1U;
 const unsigned int gGMS_MUELLER_TYPES_AVX512_MINOR = 0U;
 const unsigned int gGMS_MUELLER_TYPES_AVX512_MICRO = 0U;
 const unsigned int gGMS_MUELLER_TYPES_AVX512_FULLVER =
  1000U*gGMS_MUELLER_TYPES_AVX512_MAJOR+100U*gGMS_MUELLER_TYPES_AVX512_MINOR+10U*gGMS_MUELLER_TYPES_AVX512_MICRO;
 const char * const pgGMS_MUELLER_TYPES_AVX512_CREATION_DATE = "21-12-2021 15:03 +00200 (TUE 21 DEC 2021 15:03 GMT+2)";
 const char * const pgGMS_MUELLER_TYPES_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_MUELLER_TYPES_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_MUELLER_TYPES_AVX512_SYNOPSIS      = "AVX512 based Mueller calculus data types."


}

#include <immintrin.h>
#include "GMS_avx512c16f32.h"
#include "GMS_avx512c8f64.h"
#include "GMS_avx512vecf32.h"
#include "GMS_config.h"


namespace  gms {

         namespace math {


                     /*
                           Jones Vector based on two 16-tuple SIMD complex types.
                       */
                        typedef struct__ATTR_ALIGN__(64)  JVec2x16c16  {


                                  ZMM16c4 p;
				  ZMM16c4 s;

		       }JVec2x16c16;


		        /*
                           Jones Vector based on two 8-tuple SIMD complex types.
                       */
		       typedef struct__ATTR_ALIGN__(64)  JVec2x8c8  {
                                 
                                  ZMM8c8 p;
				  ZMM8c8 s;
		       }JVec2x8c8;


		        /*
                           Jones Matrix based on four 16-tuple SIMD complex types.
                       */
		       typedef struct  __ATTR_ALIGN__(64) JMat4x16c16 {

                                  ZMM16c4 j0;
				  ZMM16c4 j1;
				  ZMM16c4 j2;
				  ZMM16c4 j3;
		       }JMat4x16c16;

		       /*
                           Jones Matrix based on four 8-tuple SIMD complex types.
                       */
                       typedef struct __ATTR_ALIGN__(64)  JMat4x8c8 {

                                  ZMM8c8 j0;
				  ZMM8c8 j1;
				  ZMM8c8 j2;
				  ZMM8c8 j3;
		       }JMat4x8c8;


		        /*
                           Stokes Vector based on four 16-tuple SIMD real types.
                       */
		       typedef struct  __ATTR_ALIGN__(64) SVec4x16v16 {

                                  AVX512Vec16 s0;
				  AVX512Vec16 s1;
				  AVX512Vec16 s2;
				  AVX512Vec16 s3;
		       }SVec4x16v16;

		       
		         /*
                           Stokes Vector based on four 8-tuple SIMD real types.
                       */
		       typedef struct  __ATTR_ALIGN__(64) SVec4x8v8 {

                                  __m512d s0;
				  __m512d s1;
				  __m512d s2;
				  __m512d s3;
		       }SVec4x8v8;


		        /*
                           Mueller Matrix based on 16 16-tuple SIMD real types.
                       */
                       typedef struct __ATTR_ALIGN__(64)  MMat16x16v16 {

                                  AVX512Vec16 m0;
				  AVX512Vec16 m1;
				  AVX512Vec16 m2;
				  AVX512Vec16 m3;
				  AVX512Vec16 m4;
				  AVX512Vec16 m5;
				  AVX512Vec16 m6;
				  AVX512Vec16 m7;
				  AVX512Vec16 m8;
				  AVX512Vec16 m9;
				  AVX512Vec16 m10;
				  AVX512Vec16 m11;
				  AVX512Vec16 m12;
				  AVX512Vec16 m13;
				  AVX512Vec16 m14;
				  AVX512Vec16 m15;
		       }MMat16x16v16;

		       
		          /*
                           Mueller Matrix based on 16 8-tuple SIMD real types.
                       */
		       typedef struct __ATTR_ALIGN__(64)  MMat16x8v8 {

                                   __m512d m0;
				   __m512d m1;
				   __m512d m2;
				   __m512d m3;
				   __m512d m4;
				   __m512d m5;
				   __m512d m6;
				   __m512d m7;
				   __m512d m8;
				   __m512d m9;
				   __m512d m10;
				   __m512d m12;
				   __m512d m13;
				   __m512d m14;
				   __m512d m15;
		       }MMat16x8v8;

		          

		      


		      

	 
       } // math


} // gms























#endif /* __GMS_MUELLER_TYPES_AVX512_HPP__*/
