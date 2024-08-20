

#ifndef __GMS_NYSTROM_STEP_SSE_HPP__
#define __GMS_NYSTROM_STEP_SSE_HPP__ 200820240734

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


                        __ATTR_ALWAYS_INLINE__
                       	static inline
			__m128d
			nystrom_step_xmm2r8(__m128d(*f)(__m128d,
			                             __m128d),
					    __m128d y0,
					    __m128d x0,
					    const __m128d h,
					    const int32_t n) {

			    const __m128d _1_25  = _mm_set1_pd(0.04);
			    const __m128d h25    = _mm_mul_pd(_1_25,h);
			    const __m128d _1_4   = _mm_set1_pd(0.25);
			    const __m128d h4     = _mm_mul_pd(_1_4,h);
			    const __m128d _1_81  = _mm_set1_pd(0.012345679012345679012345679012);
			    const __m128d h81    = _mm_mul_pd(_1_81,h);
			    const __m128d _1_75  = _mm_set1_pd(0.013333333333333333333333333333);
			    const __m128d h75    = _mm_mul_pd(_1_45,h);
			    const __m128d _1_3   = _mm_set1_pd(0.333333333333333333333333333333);
			    const __m128d h3     = _mm_mul_pd(_1_3,h);
			    const __m128d _2_5   = _mm_set1_pd(0.4);
			    const __m128d h2_5   = _mm_mul_pd(_2_5,h);
			    const __m128d _2_3   = _mm_set1_pd(0.666666666666666666666666666667);
			    const __m128d h2_3   = _mm_mul_pd(_2_3,h);
			    const __m128d _4_5   = _mm_set1_pd(0.8);
			    const __m128d h4_5   = _mm_set1_pd(_4_5,h);
			    const __m128d _1_192 = _mm_set1_pd(0.005208333333333333333333333333);
			    const __m128d h192   = _mm_mul_pd(_1_192,h);
			    const __m128d _4     = _mm_set1_pd(4.0);
			    const __m128d _6     = _mm_set1_pd(6.0);
			    const __m128d _12    = _mm_set1_pd(12.0);
			    const __m128d _15    = _mm_set1_pd(15.0);
			    const __m128d _90    = _mm_set1_pd(90.0);
			    const __m128d _50    = _mm_set1_pd(50.0);
			    const __m128d _6     = _mm_set1_pd(6.0);
			    const __m128d _8     = _mm_set1_pd(8.0);
			    const __m128d _36    = _mm_set1_pd(36.0);
			    const __m128d _10    = _mm_set1_pd(10.0);
			    const __m128d _23    = _mm_set1_pd(23.0);
			    const __m128d _125   = _mm_set1_pd(125.0);
			    const __m128d _81    = _mm_set1_pd(81.0);
			    __m128d k1,k2,k3,k4,k5,k6;
			    __m128d t0,t1,t2,t3,t4;

			    while(--n >= 0) {
                                  k1 = f(x0,y0);
				  k2 = f(_mm_add_pd(x0,h3),
				         _mm_fmadd_pd(h3,k1,y0));
				  t0 = _mm_fmadd_pd(h25,
				                   _mm_fmadd_pd(_4,k1,
						               _mm_mul_pd(_6,k2)),y0);
				  k3 = f(_mm_add_pd(x0,h2_5),t0);
				  t1 = _mm_fmadd_pd(h4,
				                   _mm_fmadd_pd(_15,k3,
						               _mm_sub_pd(k1,
							                 _mm_mul_pd(_12,k2))),y0);
				  k4 = f(_mm_add_pd(x0,h),t1);
				  t2 = _mm_fmadd_pd(_6,k1,
				                   _mm_mul_pd(_90,k2));
				  t3 = _mm_fmadd_pd(_50,k3,
				                   _mm_mul_pd(_8,k4));
				  k5 = f(_mm_add_pd(x0,h2_3),
				         _mm_fmadd_pd(h81,_mm_sub_pd(t2,t3),y0));
				  t0 = _mm_fmadd_pd(_6,k1,
				                   _mm_fmadd_pd(_36,k2,
						               _mm_fmadd_pd(_10,k3,
							                   _mm_mul_pd(_8,k4))));
				  k6 = f(_mm_add_pd(x0,h4_5),
				         _mm_fmadd_pd(h75,t0,y0));
				  t1 = _mm_fmadd_pd(_23,k1,
				                   _mm_fmsub_pd(_125,k3,_81));
				  t2 = _mm_fmadd_pd(_125,k6,k5);
				  y0 = _mm_fmadd_pd(h192,_mm_mul_pd(t1,t2),y0);
				  x0 = _mm_add_pd(x0,h);
			    }
			    return (y0);
                      }


		        __ATTR_ALWAYS_INLINE__
                    	static inline
			__m128
			nystrom_step_xmm4r4(__m128(*f)(__m128,
			                                __m128),
					    __m128 y0,
					    __m128 x0,
					    const __m128 h,
					    const int32_t n) {

			    const __m128 _1_25  = _mm_set1_ps(0.04f);
			    const __m128 h25    = _mm_mul_ps(_1_25,h);
			    const __m128 _1_4   = _mm_set1_ps(0.25f);
			    const __m128 h4     = _mm_mul_ps(_1_4,h);
			    const __m128 _1_81  = _mm_set1_ps(0.012345679012345679012345679012f);
			    const __m128 h81    = _mm_mul_ps(_1_81,h);
			    const __m128 _1_75  = _mm_set1_ps(0.013333333333333333333333333333f);
			    const __m128 h75    = _mm_mul_ps(_1_45,h);
			    const __m128 _1_3   = _mm_set1_ps(0.333333333333333333333333333333f);
			    const __m128 h3     = _mm_mul_ps(_1_3,h);
			    const __m128 _2_5   = _mm_set1_ps(0.4f);
			    const __m128 h2_5   = _mm_mul_ps(_2_5,h);
			    const __m128 _2_3   = _mm_set1_ps(0.666666666666666666666666666667f);
			    const __m128 h2_3   = _mm_mul_ps(_2_3,h);
			    const __m128 _4_5   = _mm_set1_ps(0.8f);
			    const __m128 h4_5   = _mm_set1_ps(_4_5,h);
			    const __m128 _1_192 = _mm_set1_ps(0.005208333333333333333333333333f);
			    const __m128 h192   = _mm_mul_ps(_1_192,h);
			    const __m128 _4     = _mm_set1_ps(4.0f);
			    const __m128 _6     = _mm_set1_ps(6.0f);
			    const __m128 _12    = _mm_set1_ps(12.0f);
			    const __m128 _15    = _mm_set1_ps(15.0f);
			    const __m128 _90    = _mm_set1_ps(90.0f);
			    const __m128 _50    = _mm_set1_ps(50.0f);
			    const __m128 _6     = _mm_set1_ps(6.0f);
			    const __m128 _8     = _mm_set1_ps(8.0f);
			    const __m128 _36    = _mm_set1_ps(36.0f);
			    const __m128 _10    = _mm_set1_ps(10.0f);
			    const __m128 _23    = _mm_set1_ps(23.0f);
			    const __m128 _125   = _mm_set1_ps(125.0f);
			    const __m128 _81    = _mm_set1_ps(81.0f);
			    __m128 k1,k2,k3,k4,k5,k6;
			    __m128 t0,t1,t2,t3,t4;

			    while(--n >= 0) {
                                  k1 = f(x0,y0);
				  k2 = f(_mm_add_ps(x0,h3),
				         _mm_fmadd_ps(h3,k1,y0));
				  t0 = _mm_fmadd_ps(h25,
				                   _mm_fmadd_ps(_4,k1,
						               _mm_mul_ps(_6,k2)),y0);
				  k3 = f(_mm_add_ps(x0,h2_5),t0);
				  t1 = _mm_fmadd_ps(h4,
				                   _mm_fmadd_ps(_15,k3,
						               _mm_sub_ps(k1,
							                 _mm_mul_ps(_12,k2))),y0);
				  k4 = f(_mm_add_ps(x0,h),t1);
				  t2 = _mm_fmadd_ps(_6,k1,
				                   _mm_mul_ps(_90,k2));
				  t3 = _mm_fmadd_ps(_50,k3,
				                   _mm_mul_ps(_8,k4));
				  k5 = f(_mm_add_ps(x0,h2_3),
				         _mm_fmadd_ps(h81,_mm_sub_pd(t2,t3),y0));
				  t0 = _mm_fmadd_ps(_6,k1,
				                   _mm_fmadd_ps(_36,k2,
						               _mm_fmadd_ps(_10,k3,
							                   _mm_mul_ps(_8,k4))));
				  k6 = f(_mm_add_ps(x0,h4_5),
				         _mm_fmadd_ps(h75,t0,y0));
				  t1 = _mm_fmadd_ps(_23,k1,
				                   _mm_fmsub_ps(_125,k3,_81));
				  t2 = _mm_fmadd_ps(_125,k6,k5);
				  y0 = _mm_fmadd_ps(h192,_mm_mul_ps(t1,t2),y0);
				  x0 = _mm_add_ps(x0,h);
			    }
			    return (y0);
                      }





					    

    }


}






#endif /*__GMS_NYSTROM_STEP_SSE_HPP__*/
