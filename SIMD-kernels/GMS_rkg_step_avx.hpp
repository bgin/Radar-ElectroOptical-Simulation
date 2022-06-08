

#ifndef __GMS_RKG_STEP_AVX_HPP__
#define __GMS_RKG_STEP_AVX_HPP__ 030620221046

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
     Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_gill.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_RKG_STEP_AVX_MAJOR = 1U;
    const unsigned int GMS_RKG_STEP_AVX_MINOR = 0U;
    const unsigned int GMS_RKG_STEP_AVX_MICRO = 0U;
    const unsigned int GMS_RKG_STEP_AVX_FULLVER =
      1000U*GMS_RKG_STEP_AVX_MAJOR+
      100U*GMS_RKG_STEP_AVX_MINOR+
      10U*GMS_RKG_STEP_AVX_MICRO;
    const char * const GMS_RKG_STEP_AVX_CREATION_DATE = "03-06-2022 10:46 AM +00200 (FRI 03 JUN 2022 GMT+2)";
    const char * const GMS_RKG_STEP_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RKG_STEP_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RKG_STEP_AVX_DESCRIPTION   = "Vectorized (AVX) Runge-Kutta-Gill order 4 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"

namespace gms {

       namespace math {
////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Gill( double (*f)(double, double), double y0,          //
//                               double x0, double h, int number_of_steps );  //
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
			__m256d
			rkg_step_ymm4r8(__m256d(*f)(__m256d,
			                             __m256d),
					 __m256d y0,
					 __m256d x0,
					 const __m256d h,
					 const int32_t n) {

                           const __m256d _1_6  = _mm256_set1_pd(0.1666666666666666666667);
			   const __m256d b31   = _mm256_set1_pd(0.2071067811865475244008);
			   const __m256d b32   = _mm256_set1_pd(0.2928932188134524755992);
			   const __m256d b42   = _mm256_set1_pd(−0.7071067811865475244008);
			   const __m256d b43   = _mm256_set1_pd(1.7071067811865475244008);
			   const __m256d c2    = _mm256_set1_pd(0.5857864376269049511983);
			   const __m256d c3    = _mm256_set1_pd(3.4142135623730950488017);
			   const __m256d _1_2  = _mm256_set1_pd(0.5);
			   const __m256d h2    = _mm256_mul_pd(_1_2,h);
			   const __m256d h6    = _mm256_mul_pd(_1_6,h);
			   __m256d k1,k2,k3,k4;
			   __m256d t0,t1;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm256_add_pd(x0,h2),
				        _mm256_fmadd_pd(h2,k1,y0));
				 t0 = _mm256_fmadd_pd(h,
				                  _mm256_fmadd_pd(b31,k1,
						              _mm256_mul_pd(b32,k2)),y0);
				 k3 = f(_mm256_add_pd(x0,h2),t0);
				 x0 = _mm256_add_pd(x0,h);
				 t1 = _mm256_fmadd_pd(h,
				                  _mm256_fmadd_pd(b42,k2,
						              _mm256_mul_pd(b43,k3)),y0);
				 k4 = f(x0,t1);
				 y0 = _mm256_fmadd_pd(h6,
				                  _mm256_add_pd(_mm256_fmadd_pd(c2,k2,k1),
						                _mm256_fmadd_pd(c3,k3,k4)),y0);
			 }
			 return (y0);
		     }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256
			rkg_step_ymm8r4(__m256(*f)(__m256,
			                            __m256),
					 __m256 y0,
					 __m256 x0,
					 const __m256 h,
					 const int32_t n) {

                           const __m256 _1_6  = _mm256_set1_ps(0.1666666666666666666667F);
			   const __m256 b31   = _mm256_set1_ps(0.2071067811865475244008F);
			   const __m256 b32   = _mm256_set1_ps(0.2928932188134524755992F);
			   const __m256 b42   = _mm256_set1_ps(−0.7071067811865475244008F);
			   const __m256 b43   = _mm256_set1_ps(1.7071067811865475244008F);
			   const __m256 c2    = _mm256_set1_ps(0.5857864376269049511983F);
			   const __m256 c3    = _mm256_set1_ps(3.4142135623730950488017F);
			   const __m256 _1_2  = _mm256_set1_ps(0.5F);
			   const __m256 h2    = _mm256_mul_ps(_1_2,h);
			   const __m256 h6    = _mm256_mul_ps(_1_6,h);
			   __m256 k1,k2,k3,k4;
			   __m256 t0,t1;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm256_add_ps(x0,h2),
				        _mm256_fmadd_ps(h2,k1,y0));
				 t0 = _mm256_fmadd_ps(h,
				                  _mm256_fmadd_ps(b31,k1,
						              _mm256_mul_ps(b32,k2)),y0);
				 k3 = f(_mm256_add_ps(x0,h2),t0);
				 x0 = _mm256_add_ps(x0,h);
				 t1 = _mm256_fmadd_ps(h,
				                  _mm256_fmadd_ps(b42,k2,
						              _mm256_mul_ps(b43,k3)),y0);
				 k4 = f(x0,t1);
				 y0 = _mm256_fmadd_ps(h6,
				                  _mm256_add_ps(_mm256_fmadd_ps(c2,k2,k1),
						                _mm256_fmadd_ps(c3,k3,k4)),y0);
			 }
			 return (y0);
		     }
////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Gill_Richardson( double (*f)(double, double),          //
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


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256d
			rkg_richardson_ymm4r8(__m256d(*f)(__m256d,
			                                     __m256d),
						 __m256d y0,
						 __m256d x0,
						 const __m256d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m256d richardson[] = {_mm256_set1_pd(1.0/15.0),
			                                                    _mm256_set1_pd(1.0/31.0),
								            _mm256_set1_pd(1.0/63.0),
								            _mm256_set1_pd(1.0/127.0),
								            _mm256_set1_pd(1.0/255.0),
									    _mm256_set1_pd(1.0/511.0)};
			   constexpr int32_t MAX_COLS = 7;
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
                                     integral = rkg_step_ymm4r8(f,y0,x0,h_used,n_subints);
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
			rkg_richardson_ymm8r4(__m256(*f)(__m256,
			                                     __m256),
						 __m256 y0,
						 __m256 x0,
						 const __m256 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m256 richardson[] = { _mm256_set1_ps(1.0f/15.0f),
			                                                    _mm256_set1_ps(1.0f/31.0f),
								            _mm256_set1_ps(1.0f/63.0f),
								            _mm256_set1_ps(1.0f/127.0f),
								            _mm256_set1_ps(1.0f/255.0f),
									    _mm256_set1_ps(1.0f/511.0f)};
			   constexpr int32_t MAX_COLS = 7;
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
                                     integral = rkg_step_ymm8r4(f,y0,x0,h_used,n_subints);
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



   }


}











#endif /*__GMS_RKG_STEP_AVX_HPP__*/
