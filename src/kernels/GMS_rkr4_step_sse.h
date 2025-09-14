

#ifndef __GMS_RKR4_STEP_SSE_H__
#define __GMS_RKR4_STEP_SSE_H__ 050920240819

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
     Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_ralston_4.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

#if !defined(__FMA__) || !defined(__AVX512F__)
#error "Support of FMA ISA is required!!"
#endif


namespace file_version {

    const unsigned int GMS_RKR4_STEP_SSE_MAJOR = 1U;
    const unsigned int GMS_RKR4_STEP_SSE_MINOR = 0U;
    const unsigned int GMS_RKR4_STEP_SSE_MICRO = 0U;
    const unsigned int GMS_RKR4_STEP_SSE_FULLVER =
      1000U*GMS_RKR4_STEP_SSE_MAJOR+
      100U*GMS_RKR4_STEP_SSE_MINOR+
      10U*GMS_RKR4_STEP_SSE_MICRO;
    const char * const GMS_RKR4_STEP_SSE_CREATION_DATE = "05-09-2024 08:19 PM +00200 (THR 05 SEP 2024 GMT+2)";
    const char * const GMS_RKR4_STEP_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RKR4_STEP_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RKR4_STEP_SSE_DESCRIPTION   = "Vectorized (SSE) Runge-Kutta-Ralston order 4 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {


           namespace math {

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The fourth order Ralston's method is a four stage Runge-Kutta method   //
//     for approximating the solution of the differential equation            //
//     y'(x) = f(x,y) with initial condition y(x0) = c.  The fourth order     //
//     Ralston's method evaluates f(x,y) four times per step. For step i+1,   //
//     y[i+1] = y[i] + (263 + 24s5)/1812  k1 + (125-1000s5)/3828 k2           //
//                   + 1024*(3346 + 1623s5)/5924787 k3 + (30 - 4s5)/123 k4,   //
//     where                                                                  //
//     k1 = h * f(x[i],y[i]),                                                 //
//     k2 = h * f(x[i]+2/5 h, y[i]+2/5 k1),                                   //
//     k3 = h * f(x[i] + (7/8 - 3s5/16)h, y[i] + (-2889+1428s5)/1024 k1       //
//                                             + (3785-1620s5)/1024 k2)       //
//     k4 = h * f(x[i] +h, y[i] + (-3365+2094s5)/6040 k1                      //
//                                (-975-3046s5)/2552 k2                       //
//                                (467040+203968s5)/2400845 k3) )             //
//     and s5 = sqrt(5), x[i] = x0 + i * h.                                   //
//                                                                            //
//     This version of the Runge-Kutta method is a fourth order method.       //
//     Richardson extrapolation can be used to increase the order and         //
//     accuracy.                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


                        __ATTR_HOT__
                        __ATTR_VECTORCALL__
			__m128d
			rkr4_step_xmm2r8(__m128d(*f)(__m128d,
			                             __m128d),
					 __m128d y0,
					 __m128d x0,
					 const __m128d h,
					 const int32_t n); 


		         __ATTR_HOT__
                        __ATTR_VECTORCALL__
			__m128
			rkr4_step_xmm4r4(__m128(*f)(__m128,
			                             __m128),
					 __m128 y0,
					 __m128 x0,
					 const __m128 h,
					 const int32_t n); 

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Ralston_4_Richardson( double (*f)(double, double),     //
//                       double y0, double x0, double h, int number_of_steps, //
//                                                   int richardson_columns)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the 4th order Runge-Kutta method described above     //
//     together with Richardson extrapolation to approximate the solution at  //
//     x = x0 + h * number_of_steps of the initial value problem y'=f(x,y),   //
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
//            The number of steps. Must be nonnegative.                       //
//     int    richardson_columns                                              //
//            The maximum number of columns to use in the Richardson          //
//            extrapolation to the limit.                                     //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <algorithm>


                        __ATTR_HOT__
                        __ATTR_VECTORCALL__
			__m128d
			rkr4_richardson_xmm2r8(__m128d(*f)(__m128d,
			                                     __m128d),
						 __m128d y0,
						 __m128d x0,
						 const __m128d h,
						 const int32_t n,
						 int32_t n_cols); 


		        __ATTR_HOT__
                        __ATTR_VECTORCALL__
			__m128
			rkr4_richardson_xmm4r4(__m128(*f)(__m128,
			                                     __m128),
						 __m128 y0,
						 __m128 x0,
						 const __m128 h,
						 const int32_t n,
						 int32_t n_cols); 
	   
     }


}






#endif /*__GMS_RKR4_STEP_SSE_H__*/
