

#ifndef __GMS_NYSTROM_STEP_SSE_H__
#define __GMS_NYSTROM_STEP_SSE_H__ 200820240734

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
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_nystrom.html
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_NYSTROM_STEP_SSE_MAJOR = 1U;
    const unsigned int GMS_NYSTROM_STEP_SSE_MINOR = 0U;
    const unsigned int GMS_NYSTROM_STEP_SSE_MICRO = 0U;
    const unsigned int GMS_NYSTROM_STEP_SSE_FULLVER =
      1000U*GMS_NYSTROM_STEP_SSE_MAJOR+
      100U*GMS_NYSTROM_STEP_SSE_MINOR+
      10U*GMS_NYSTROM_STEP_SSE_MICRO;
    const char * const GMS_NYSTROM_STEP_SSE_CREATION_DATE = "20-08-2024 07:34 PM +00200 (TUE 20 AUG 2024 GMT+2)";
    const char * const GMS_NYSTROM_STEPSSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_NYSTROM_STEP_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_NYSTROM_STEP_SSE_DESCRIPTION   = "Vectorized (SSE) Runge-Kutta-Nystrom order 5 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {


      namespace math {

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Nystrom( double (*f)(double, double), double y0,       //
//                               double x0, double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Nystrom method described above to    //
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


                     
			__m128d
			nystrom_step_xmm2r8(__m128d(*f)(__m128d,
			                             __m128d),
					    __m128d y0,
					    __m128d x0,
					    const __m128d h,
					    const int32_t n); 
					    

		        
			__m128
			nystrom_step_xmm4r4(__m128(*f)(__m128,
			                                __m128),
					    __m128 y0,
					    __m128 x0,
					    const __m128 h,
					    const int32_t n); 




					    

    }


}






#endif /*__GMS_NYSTROM_STEP_SSE_H__*/
