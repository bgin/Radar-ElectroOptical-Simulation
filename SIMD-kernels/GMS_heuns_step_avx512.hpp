

#ifndef __GMS_HEUNS_STEP_AVX512_HPP__
#define __GMS_HEUNS_STEP_AVX512_HPP__ 010620221256

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

/*
    Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/heuns_method.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_HEUNS_STEP_AVX512_MAJOR = 1U;
    const unsigned int GMS_HEUNS_STEP_AVX512_MINOR = 0U;
    const unsigned int GMS_HEUNS_STEP_AVX512_MICRO = 0U;
    const unsigned int GMS_HEUNS_STEP_AVX512_FULLVER =
      1000U*GMS_HEUNS_STEP_AVX512_MAJOR+
      100U*GMS_HEUNS_STEP_AVX512_MINOR+
      10U*GMS_HEUNS_STEP_AVX512_MICRO;
    const char * const GMS_HEUNS_STEP_AVX512_CREATION_DATE = "01-06-2022 12:57 PM +00200 (WED 01 JUN 2022 GMT+2)";
    const char * const GMS_HEUNS_STEPAVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_HEUNS_STEP_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_HEUNS_STEP_AVX512_DESCRIPTION   = "Vectorized (AVX512) Heuns step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {


         namespace  math {

////////////////////////////////////////////////////////////////////////////////
//  double Heuns_Method( double (*f)(double, double), double y0, double x0,   //
//                                          double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Heun's method to approximate the solution at x =     //
//     x0 + h * number_of_steps of the initial value problem y'=f(x,y),       //
//     y(x0) = y0.                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     double y0                                                              //
//            The initial value of y at x = x0.                               //
//     double x0                                                              //
//            The initial value of x.                                         //
//     double h                                                               //
//            The step size.                                                  //
//     int    number_of_steps                                                 //
//            The number of steps. Must be a nonnegative integer.             //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//

                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d
			heuns_step_zmm8r8(__m512d(*f)(__m512d,
			                              __m512d),
					  __m512d y0,
					  __m512d x0,
					  __m512d h,
					  const int32_t n) {

                          const __m512d _0_5 = _mm512_set1_pd(0.5);
			  __m512d k1;
			  while(--n >= 0) {
                                k1 = _mm512_mul_pd(h,f(x0,y0));
				x0 = _mm512_add_pd(x0,h);
				y0 = _mm512_add_pd(y0,_mm512_mul_pd(_0_5,
				                                _mm512_fmadd_pd(h,f(x0,_mm512_add_pd(y0,1),k1))));
			  }
			  return (y0);
		      }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			heuns_step_zmm16r4(__m512(*f)(__m512,
			                              __m512),
					  __m512 y0,
					  __m512 x0,
					  __m512 h,
					  const int32_t n) {

                          const __m512 _0_5 = _mm512_set1_ps(0.5f);
			  __m512 k1;
			  while(--n >= 0) {
                                k1 = _mm512_mul_ps(h,f(x0,y0));
				x0 = _mm512_add_ps(x0,h);
				y0 = _mm512_add_ps(y0,_mm512_mul_ps(_0_5,
				                                _mm512_fmadd_ps(h,f(x0,_mm512_add_ps(y0,1),k1))));
			  }
			  return (y0);
		      }


      } //math

 
} //gms




















#endif /*__GMS_HEUNS_STEP_AVX512_HPP__*/
