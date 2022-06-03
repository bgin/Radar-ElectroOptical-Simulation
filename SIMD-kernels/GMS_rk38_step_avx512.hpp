

#ifndef __GMS_RK38_STEP_AVX512_HPP__
#define __GMS_RK38_STEP_AVX512_HPP__  030620220818

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
     Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_3_8.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_RK38_STEP_AVX512_MAJOR = 1U;
    const unsigned int GMS_RK38_STEP_AVX512_MINOR = 0U;
    const unsigned int GMS_RK38_STEP_AVX512_MICRO = 0U;
    const unsigned int GMS_RK38_STEP_AVX512_FULLVER =
      1000U*GMS_RK38_STEP_AVX512_MAJOR+
      100U*GMS_RK38_STEP_AVX512_MINOR+
      10U*GMS_RK38_STEP_AVX512_MICRO;
    const char * const GMS_RK38_STEP_AVX512_CREATION_DATE = "03-06-2022 08:18 PM +00200 (FRI 03 JUN 2022 GMT+2)";
    const char * const GMS_RK38_STEPAVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RK38_STEP_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RK38_STEP_AVX512_DESCRIPTION   = "Vectorized (AVX512) Runge-Kutta-Ralston order 3/8 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {


     namespace math {

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_3_8( double (*f)(double, double), double y0, double x0,//
//                                          double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the 4th order Runge-Kutta method described above to  //
//     approximate the solution at x = x0 + h * number_of_steps of the initial//
//     value problem y'=f(x,y), y(x0) = y0.                                   //
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


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d
			rk38_step_zmm8r8(__m512d(*f)(__m512d,
			                             __m512d),
					 __m512d y0,
					 __m512d x0,
					 const __m512d h,
					 const int32_t n) {

                           const __m512d _1_8   = _mm512_set1_pd(0.125);
			   const __m512d _1_3   = _mm512_set1_pd(0.3333333333333333333333333);
			   const __m512d _2_3   = _mm512_set1_pd(0.6666666666666666666667);
			   const __m512d _3     = _mm512_set1_pd(3.0);
			   const __m512d h13    = _mm512_mul_pd(_1_3,h);
			   const __m512d h23    = _mm512_mul_pd(_2_3,h);
			   const __m512d h18    = _mm512_mul_pd(_1_8,h); 
			   __m512d k1,k2,k3,k4;
			   __m512d t0,t1;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm512_add_pd(x0,h13),
				        _mm512_fmadd_pd(h13,k1,y0));
				 t0 = _mm512_fmadd_pd(h,
				                  _mm512_fmsub_pd(_1_3,k1,k2),y0);
				 k3 = f(_mm512_add_pd(x0,h23),t0);
				 x0 = _mm512_add_pd(x0,h);
				 t1 = _mm512_fmadd_pd(h,
				                  _mm512_sub_pd(k1,
						            _mm512_add_pd(k2,k3)),y0);
				 k4 = f(x0,t1);
				 y0 = _mm512_fmadd_pd(h18,
				                  _mm512_add_pd(
						            _mm512_fmadd_pd(_3,k2,k1),
							    _mm512_fmadd_pd(_3,k3,k4)),y0);
			   }
			   return (y0);
		    }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			rk38_step_zmm16r4(__m512(*f)(__m512,
			                             __m512),
					 __m512 y0,
					 __m512 x0,
					 const __m512 h,
					 const int32_t n) {

                           const __m512 _1_8   = _mm512_set1_ps(0.125f);
			   const __m512 _1_3   = _mm512_set1_ps(0.3333333333333333333333333f);
			   const __m512 _2_3   = _mm512_set1_ps(0.6666666666666666666667f);
			   const __m512 _3     = _mm512_set1_ps(3.0f);
			   const __m512 h13    = _mm512_mul_ps(_1_3,h);
			   const __m512 h23    = _mm512_mul_ps(_2_3,h);
			   const __m512 h18    = _mm512_mul_ps(_1_8,h); 
			   __m512 k1,k2,k3,k4;
			   __m512 t0,t1;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm512_add_ps(x0,h13),
				        _mm512_fmadd_ps(h13,k1,y0));
				 t0 = _mm512_fmadd_ps(h,
				                  _mm512_fmsub_ps(_1_3,k1,k2),y0);
				 k3 = f(_mm512_add_ps(x0,h23),t0);
				 x0 = _mm512_add_ps(x0,h);
				 t1 = _mm512_fmadd_ps(h,
				                  _mm512_sub_ps(k1,
						            _mm512_add_ps(k2,k3)),y0);
				 k4 = f(x0,t1);
				 y0 = _mm512_fmadd_ps(h18,
				                  _mm512_add_ps(
						            _mm512_fmadd_ps(_3,k2,k1),
							    _mm512_fmadd_ps(_3,k3,k4)),y0);
			   }
			   return (y0);
		    }



		    

   }



}









#endif /*__GMS_RK38_STEP_AVX512_HPP__*/
