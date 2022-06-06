

#ifndef __GMS_EULER_STEP_AVX_HPP__
#define __GMS_EULER_STEP_AVX_HPP__ 010620221448

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
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/improved_euler_method.html
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_EULER_STEP_AVX_MAJOR = 1U;
    const unsigned int GMS_EULER_STEP_AVX_MINOR = 0U;
    const unsigned int GMS_EULER_STEP_AVX_MICRO = 0U;
    const unsigned int GMS_EULER_STEP_AVX_FULLVER =
      1000U*GMS_EULER_STEP_AVX_MAJOR+
      100U*GMS_EULER_STEP_AVX_MINOR+
      10U*GMS_EULER_STEP_AVX_MICRO;
    const char * const GMS_EULER_STEP_AVX_CREATION_DATE = "01-06-2022 14:48 PM +00200 (WED 01 JUN 2022 GMT+2)";
    const char * const GMS_EULER_STEP_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EULER_STEP_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EULER_STEP_AVX_DESCRIPTION   = "Vectorized (AVX) Runge-Kutta order 4 step."

}


#include <immintrin.h>
#include "GMS_config.h"

namespace  gms {

          namespace math {

////////////////////////////////////////////////////////////////////////////////
//  double Improved_Euler_Method( double (*f)(double, double), double y0,     //
//                               double x0, double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the improved Euler method to approximate the solution//
//     at x = x0 + h * number_of_steps of the initial value problem y'=f(x,y) //
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


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256d
			euler_step_ymm4r8(__m256d(*f)(__m256d,
			                              __m256d),
					  __m256d y0,
					  __m256d x0,
					  const __m256d h,
					  const int32_t n) {

                           const __m256d _0_5 = _mm256_set1_pd(0.5);
			   const __m256d h2   = _mm256_mul_pd(_0_5,h);
			   __m512d k1;
			   while(--n >= 0) {
                                 k1 = _mm256_mul_pd(h,f(x0,y0));
				 y0 = _mm256_add_pd(y0,_mm256_mul_pd(h,
				                                 f(_mm256_add_pd(x0,h2),
								   _mm256_add_pd(y0,k1))));
				 x0 = _mm256_add_pd(x0,h);
			   }
			   return (y0);
		      }


		      
	               __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256
			euler_step_ymm8r4(__m256(*f)(__m256,
			                             __m256),
					  __m256 y0,
					  __m256 x0,
					  const __m256 h,
					  const int32_t n) {

                           const __m256 _0_5 = _mm256_set1_ps(0.5f);
			   const __m256 h2   = _mm256_mul_ps(_0_5,h);
			   __m512 k1;
			   while(--n >= 0) {
                                 k1 = _mm256_mul_ps(h,f(x0,y0));
				 y0 = _mm256_add_ps(y0,_mm256_mul_ps(h,
				                                 f(_mm256_add_ps(x0,h2),
								   _mm256_add_ps(y0,k1))));
				 x0 = _mm256_add_ps(x0,h);
			   }
			   return (y0);
		      }


	    ////////////////////////////////////////////////////////////////////////////////
//  double Improved_Euler_Method_Richardson( double (*f)(double, double),     //
//                       double y0, double x0, double h, int number_of_steps, //
//                                                   int richardson_columns)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the improved Euler method together with Richardson   //
//     extrapolation to approximate the solution at x = x0 + h*number_of_steps//
//     of the initial value problem y'=f(x,y), y(x0) = y0.                    //
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


                  
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256d
			euler_richardson_ymm4r8(__m512d(*f)(__m256d,
			                                     __m256d),
						 __m256d y0,
						 __m256d x0,
						 const __m256d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m256d richardson[] = {_mm256_set1_pd(1.0/3.0),
			                                                    _mm256_set1_pd(1.0/7.0),
								            _mm256_set1_pd(1.0/15.0),
								            _mm256_set1_pd(1.0/31.0),
								            _mm256_set1_pd(1.0/63.0)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(32) __m256d dt[MAX_COLS];
			   const __m256d _1_2 = _mm256_set1_pd(0.5);
			   __m256d integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = euler_step_ymm4r8(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm256_sub_pd(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm256_fmadd_pd(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm256_mul_pd(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm256_add_pd(x0,h);
			   }
			   return (y0);
		    }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256
			euler_richardson_ymm8r4(__m256(*f)(__m256,
			                                     __m256),
						 __m256 y0,
						 __m256 x0,
						 const __m256 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m256 richardson[] = { _mm256_set1_ps(1.0f/3.0f),
			                                                    _mm256_set1_ps(1.0f/7.0f),
								            _mm256_set1_ps(1.0f/15.0f),
								            _mm256_set1_ps(1.0f/31.0f),
								            _mm256_set1_ps(1.0f/63.0f)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(32) __m256 dt[MAX_COLS];
			   const __m256 _1_2 = _mm256_set1_ps(0.5f);
			   __m256 integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = euler_step_ymm8r4(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm256_sub_ps(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm256_fmadd_ps(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm256_mul_ps(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm256_add_ps(x0,h);
			   }
			   return (y0);
		    }



    } // math

} //gms





#endif /*__GMS_EULER_STEP_AVX_HPP__*/
