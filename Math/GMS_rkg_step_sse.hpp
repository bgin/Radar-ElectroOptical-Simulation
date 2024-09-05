

#ifndef __GMS_RKG_STEP_SSE_HPP__
#define __GMS_RKG_STEP_SSE_HPP__ 030620221046

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

#if !defined(__FMA__) || !defined(__AVX512F__)
#error "Support of FMA ISA is required!!"
#endif

namespace file_version {

    const unsigned int GMS_RKG_STEP_SSE_MAJOR = 1U;
    const unsigned int GMS_RKG_STEP_SSE_MINOR = 0U;
    const unsigned int GMS_RKG_STEP_SSE_MICRO = 0U;
    const unsigned int GMS_RKG_STEP_SSE_FULLVER =
      1000U*GMS_RKG_STEP_SSE_MAJOR+
      100U*GMS_RKG_STEP_SSE_MINOR+
      10U*GMS_RKG_STEP_SSE_MICRO;
    const char * const GMS_RKG_STEP_SSE_CREATION_DATE = "03-06-2022 10:46 AM +00200 (FRI 03 JUN 2022 GMT+2)";
    const char * const GMS_RKG_STEP_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RKG_STEP_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RKG_STEP_SSE_DESCRIPTION   = "Vectorized (SSE) Runge-Kutta-Gill order 4 step."

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
                       
	                static inline
			__m128d
			rkg_step_xmm2r8(__m128d(*f)(__m128d,
			                             __m128d),
					 __m128d y0,
					 __m128d x0,
					 const __m128d h,
					 const int32_t n) {

                           const __m128d _1_6  = _mm_set1_pd(0.1666666666666666666667);
			   const __m128d b31   = _mm_set1_pd(0.2071067811865475244008);
			   const __m128d b32   = _mm_set1_pd(0.2928932188134524755992);
			   const __m128d b42   = _mm_set1_pd(−0.7071067811865475244008);
			   const __m128d b43   = _mm_set1_pd(1.7071067811865475244008);
			   const __m128d c2    = _mm_set1_pd(0.5857864376269049511983);
			   const __m128d c3    = _mm_set1_pd(3.4142135623730950488017);
			   const __m128d _1_2  = _mm_set1_pd(0.5);
			   const __m128d h2    = _mm_mul_pd(_1_2,h);
			   const __m128d h6    = _mm_mul_pd(_1_6,h);
			   __m128d k1,k2,k3,k4;
			   __m128d t0,t1;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm_add_pd(x0,h2),
				        _mm_fmadd_pd(h2,k1,y0));
				 t0 = _mm_fmadd_pd(h,
				                  _mm_fmadd_pd(b31,k1,
						              _mm_mul_pd(b32,k2)),y0);
				 k3 = f(_mm_add_pd(x0,h2),t0);
				 x0 = _mm_add_pd(x0,h);
				 t1 = _mm_fmadd_pd(h,
				                  _mm_fmadd_pd(b42,k2,
						              _mm_mul_pd(b43,k3)),y0);
				 k4 = f(x0,t1);
				 y0 = _mm_fmadd_pd(h6,
				                  _mm_add_pd(_mm_fmadd_pd(c2,k2,k1),
						                _mm_fmadd_pd(c3,k3,k4)),y0);
			 }
			 return (y0);
		     }


		        __ATTR_ALWAYS_INLINE__
                       	static inline
			__m128
			rkg_step_xmm4r4(__m128(*f)(__m128,
			                            __m128),
					 __m128 y0,
					 __m128 x0,
					 const __m128 h,
					 const int32_t n) {

                           const __m128 _1_6  = _mm_set1_ps(0.1666666666666666666667F);
			   const __m128 b31   = _mm_set1_ps(0.2071067811865475244008F);
			   const __m128 b32   = _mm_set1_ps(0.2928932188134524755992F);
			   const __m128 b42   = _mm_set1_ps(−0.7071067811865475244008F);
			   const __m128 b43   = _mm_set1_ps(1.7071067811865475244008F);
			   const __m128 c2    = _mm_set1_ps(0.5857864376269049511983F);
			   const __m128 c3    = _mm_set1_ps(3.4142135623730950488017F);
			   const __m128 _1_2  = _mm_set1_ps(0.5F);
			   const __m128 h2    = _mm_mul_ps(_1_2,h);
			   const __m128 h6    = _mm_mul_ps(_1_6,h);
			   __m128 k1,k2,k3,k4;
			   __m128 t0,t1;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm_add_ps(x0,h2),
				        _mm_fmadd_ps(h2,k1,y0));
				 t0 = _mm_fmadd_ps(h,
				                  _mm_fmadd_ps(b31,k1,
						              _mm_mul_ps(b32,k2)),y0);
				 k3 = f(_mm_add_ps(x0,h2),t0);
				 x0 = _mm_add_ps(x0,h);
				 t1 = _mm_fmadd_ps(h,
				                  _mm_fmadd_ps(b42,k2,
						              _mm_mul_ps(b43,k3)),y0);
				 k4 = f(x0,t1);
				 y0 = _mm_fmadd_ps(h6,
				                  _mm_add_ps(_mm_fmadd_ps(c2,k2,k1),
						                _mm_fmadd_ps(c3,k3,k4)),y0);
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
                       	static inline
			__m128d
			rkg_richardson_xmm2r8(__m128d(*f)(__m128d,
			                                     __m128d),
						 __m128d y0,
						 __m128d x0,
						 const __m128d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(16) const __m128d richardson[] = {_mm_set1_pd(1.0/15.0),
			                                                    _mm_set1_pd(1.0/31.0),
								            _mm_set1_pd(1.0/63.0),
								            _mm_set1_pd(1.0/127.0),
								            _mm_set1_pd(1.0/255.0),
									    _mm_set1_pd(1.0/511.0)};
			   constexpr int32_t MAX_COLS = 7;
			   __ATTR_ALIGN__(16) __m128d dt[MAX_COLS];
			   const __m128d _1_2 = _mm_set1_pd(0.5);
			   __m128d integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = rkg_step_xmm2r8(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm_sub_pd(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm_fmadd_pd(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm_mul_pd(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm_add_pd(x0,h);
			   }
			   return (y0);
		    }


		        __ATTR_ALWAYS_INLINE__
                       	static inline
			__m128
			rkg_richardson_xmm4r4(__m128(*f)(__m128,
			                                     __m128),
						 __m128 y0,
						 __m128 x0,
						 const __m128 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(16) const __m128 richardson[] = { _mm_set1_ps(1.0f/15.0f),
			                                                    _mm_set1_ps(1.0f/31.0f),
								            _mm_set1_ps(1.0f/63.0f),
								            _mm_set1_ps(1.0f/127.0f),
								            _mm_set1_ps(1.0f/255.0f),
									    _mm_set1_ps(1.0f/511.0f)};
			   constexpr int32_t MAX_COLS = 7;
			   __ATTR_ALIGN__(16) __m128 dt[MAX_COLS];
			   const __m128 _1_2 = _mm_set1_ps(0.5f);
			   __m128 integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = rkg_step_xmm4r4(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm_sub_ps(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm_fmadd_ps(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm_mul_ps(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm_add_ps(x0,h);
			   }
			   return (y0);
		    }



   }


}











#endif /*__GMS_RKG_STEP_SSE_HPP__*/
