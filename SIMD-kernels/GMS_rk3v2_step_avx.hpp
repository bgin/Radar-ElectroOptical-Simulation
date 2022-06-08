

#ifndef __GMS_RK3V2_STEP_AVX_HPP__
#define __GMS_RK3V2_STEP_AVX_HPP__ 020620221301



/*
     Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_v2_3.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_RK3V2_STEP_AVX_MAJOR = 1U;
    const unsigned int GMS_RK3V2_STEP_AVX_MINOR = 0U;
    const unsigned int GMS_RK3V2_STEP_AVX_MICRO = 0U;
    const unsigned int GMS_RK3V2_STEP_AVX_FULLVER =
      1000U*GMS_RK3V2_STEP_AVX_MAJOR+
      100U*GMS_RK3V2_STEP_AVX_MINOR+
      10U*GMS_RK3V2_STEP_AVX_MICRO;
    const char * const GMS_RK3V2_STEP_AVX_CREATION_DATE = "02-06-2022 13:01 PM +00200 (THR 02 JUN 2022 GMT+2)";
    const char * const GMS_RK3V2_STEP_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RK3V2_STEP_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RK3V2_STEP_AVX_DESCRIPTION   = "Vectorized (AVX) Runge-Kutta order 3 (version 2) step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"



namespace  gms {


            namespace math {

	    
////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_v2_3( double (*f)(double, double), double y0,          //
//                               double x0, double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the third order Runge-Kutta method described above   //
//     to approximate the solution at x = x0 + h * number_of_steps of the     //
//     initial value problem y'=f(x,y), y(x0) = y0.                           //
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
			rk3v2_step_ymm4r8(__m256d(*f)(__m256d,
			                             __m256d),
					  __m256d y0,
					  __m256d x0,
					  const __m256d h,
					  const int32_t n) {

                           const __m256d _1_3 = _mm256_set1_pd(0.3333333333333333333333333333333);
			   const __m256d _2_3 = _mm256_set1_pd(0.6666666666666666666667);
			   const __m256d _1_4 = _mm256_set1_pd(0.25);
			   const __m256d _3   = _mm256_set1_pd(3.0);
			   const __m256d h13  = _mm256_mul_pd(_1_3,h);
			   const __m256d h23  = _mm256_mul_pd(_2_3,h);
			   const __m256d h4   = _mm256_mul_pd(_1_4,h);
			   __m256d k1,k2;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm256_add_pd(x0,h13),
				        _mm256_fmadd_pd(h13,k1,y0));
				 y0 = _mm256_add_pd(y0,_mm256_mul_pd(h4,
				                                _mm256_fmadd_pd(_3,f(_mm256_add_pd(x0,h23),
								                     _mm256_fmadd_pd(h23,k2,y0)),k1)));
				 x0 = _mm256_add_pd(x0,h);
			   }
			   return (y0);
		      }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			rk3v2_step_ymm8r4(__m512(*f)(__m512,
			                             __m512),
					  __m512 y0,
					  __m512 x0,
					  const __m512 h,
					  const int32_t n) {

                           const __m512 _1_3 = _mm256_set1_ps(0.3333333333333333333333333333333f);
			   const __m512 _2_3 = _mm256_set1_ps(0.6666666666666666666667f);
			   const __m512 _1_4 = _mm256_set1_ps(0.25f);
			   const __m512 _3   = _mm256_set1_ps(3.0f);
			   const __m512 h13  = _mm256_mul_ps(_1_3,h);
			   const __m512 h23  = _mm256_mul_ps(_2_3,h);
			   const __m512 h4   = _mm256_mul_ps(_1_4,h);
			   __m512 k1,k2;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm256_add_ps(x0,h13),
				        _mm256_fmadd_ps(h13,k1,y0));
				 y0 = _mm256_add_ps(y0,_mm256_mul_ps(h4,
				                                _mm256_fmadd_ps(_3,f(_mm256_add_ps(x0,h23),
								                     _mm256_fmadd_ps(h23,k2,y0)),k1)));
				 x0 = _mm256_add_ps(x0,h);
			   }
			   return (y0);
		      }

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_v2_3_Richardson( double (*f)(double, double),          //
//                       double y0, double x0, double h, int number_of_steps, //
//                                                   int richardson_columns)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the third order Runge-Kutta method described above   //
//     together with Richardson extrapolation to approximate the solution     //
//     at x = x0 + h * number_of_steps of the initial value problem           //
//     y'=f(x,y), y(x0) = y0.                                                 //
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
			rk3v2_richardson_ymm4r8(__m256d(*f)(__m256d,
			                                     __m256d),
						 __m256d y0,
						 __m256d x0,
						 const __m256d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(64) const __m256d richardson[] = {_mm256_set1_pd(1.0/7.0),
			                                                    _mm256_set1_pd(1.0/15.0),
								            _mm256_set1_pd(1.0/31.0),
								            _mm256_set1_pd(1.0/63.0),
								            _mm256_set1_pd(1.0/127.0)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(64) __m256d dt[MAX_COLS];
			   const __m256d _1_2 = _mm256_set1_pd(0.5);
			   __m256d integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = rk3v2_step_ymm4r8(f,y0,x0,h_used,n_subints);
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
			__m512
			rk3v2_richardson_ymm8r4(__m512(*f)(__m512,
			                                     __m512),
						 __m512 y0,
						 __m512 x0,
						 const __m512 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(64) const __m512 richardson[] = { _mm256_set1_ps(1.0f/7.0f),
			                                                    _mm256_set1_ps(1.0f/15.0f),
								            _mm256_set1_ps(1.0f/31.0f),
								            _mm256_set1_ps(1.0f/63.0f),
								            _mm256_set1_ps(1.0f/127.0f)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(64) __m512 dt[MAX_COLS];
			   const __m512 _1_2 = _mm256_set1_ps(0.5f);
			   __m512 integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = rk3v2_step_ymm8r4(f,y0,x0,h_used,n_subints);
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







#endif /*__GMS_RK3V2_STEP_AVX_HPP__*/
