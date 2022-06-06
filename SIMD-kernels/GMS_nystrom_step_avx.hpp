

#ifndef __GMS_NYSTROM_STEP_AVX_HPP__
#define __GMS_NYSTROM_STEP_AVX_HPP__ 040620220959

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

    const unsigned int GMS_NYSTROM_STEP_AVX_MAJOR = 1U;
    const unsigned int GMS_NYSTROM_STEP_AVX_MINOR = 0U;
    const unsigned int GMS_NYSTROM_STEP_AVX_MICRO = 0U;
    const unsigned int GMS_NYSTROM_STEP_AVX_FULLVER =
      1000U*GMS_NYSTROM_STEP_AVX_MAJOR+
      100U*GMS_NYSTROM_STEP_AVX_MINOR+
      10U*GMS_NYSTROM_STEP_AVX_MICRO;
    const char * const GMS_NYSTROM_STEP_AVX_CREATION_DATE = "04-06-2022 09:59 PM +00200 (SAT 04 JUN 2022 GMT+2)";
    const char * const GMS_NYSTROM_STEPAVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_NYSTROM_STEP_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_NYSTROM_STEP_AVX_DESCRIPTION   = "Vectorized (AVX) Runge-Kutta-Nystrom order 5 step."

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
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m256d
			nystrom_step_ymm4r8(__m256d(*f)(__m256d,
			                             __m256d),
					    __m256d y0,
					    __m256d x0,
					    const __m256d h,
					    const int32_t n) {

			    const __m256d _1_25  = _mm256_set1_pd(0.04);
			    const __m256d h25    = _mm256_mul_pd(_1_25,h);
			    const __m256d _1_4   = _mm256_set1_pd(0.25);
			    const __m256d h4     = _mm256_mul_pd(_1_4,h);
			    const __m256d _1_81  = _mm256_set1_pd(0.012345679012345679012345679012);
			    const __m256d h81    = _mm256_mul_pd(_1_81,h);
			    const __m256d _1_75  = _mm256_set1_pd(0.013333333333333333333333333333);
			    const __m256d h75    = _mm256_mul_pd(_1_45,h);
			    const __m256d _1_3   = _mm256_set1_pd(0.333333333333333333333333333333);
			    const __m256d h3     = _mm256_mul_pd(_1_3,h);
			    const __m256d _2_5   = _mm256_set1_pd(0.4);
			    const __m256d h2_5   = _mm256_mul_pd(_2_5,h);
			    const __m256d _2_3   = _mm256_set1_pd(0.666666666666666666666666666667);
			    const __m256d h2_3   = _mm256_mul_pd(_2_3,h);
			    const __m256d _4_5   = _mm256_set1_pd(0.8);
			    const __m256d h4_5   = _mm256_set1_pd(_4_5,h);
			    const __m256d _1_192 = _mm256_set1_pd(0.005208333333333333333333333333);
			    const __m256d h192   = _mm256_mul_pd(_1_192,h);
			    const __m256d _4     = _mm256_set1_pd(4.0);
			    const __m256d _6     = _mm256_set1_pd(6.0);
			    const __m256d _12    = _mm256_set1_pd(12.0);
			    const __m256d _15    = _mm256_set1_pd(15.0);
			    const __m256d _90    = _mm256_set1_pd(90.0);
			    const __m256d _50    = _mm256_set1_pd(50.0);
			    const __m256d _6     = _mm256_set1_pd(6.0);
			    const __m256d _8     = _mm256_set1_pd(8.0);
			    const __m256d _36    = _mm256_set1_pd(36.0);
			    const __m256d _10    = _mm256_set1_pd(10.0);
			    const __m256d _23    = _mm256_set1_pd(23.0);
			    const __m256d _125   = _mm256_set1_pd(125.0);
			    const __m256d _81    = _mm256_set1_pd(81.0);
			    __m256d k1,k2,k3,k4,k5,k6;
			    __m256d t0,t1,t2,t3,t4;

			    while(--n >= 0) {
                                  k1 = f(x0,y0);
				  k2 = f(_mm256_add_pd(x0,h3),
				         _mm256_fmadd_pd(h3,k1,y0));
				  t0 = _mm256_fmadd_pd(h25,
				                   _mm256_fmadd_pd(_4,k1,
						               _mm256_mul_pd(_6,k2)),y0);
				  k3 = f(_mm256_add_pd(x0,h2_5),t0);
				  t1 = _mm256_fmadd_pd(h4,
				                   _mm256_fmadd_pd(_15,k3,
						               _mm256_sub_pd(k1,
							                 _mm256_mul_pd(_12,k2))),y0);
				  k4 = f(_mm256_add_pd(x0,h),t1);
				  t2 = _mm256_fmadd_pd(_6,k1,
				                   _mm256_mul_pd(_90,k2));
				  t3 = _mm256_fmadd_pd(_50,k3,
				                   _mm256_mul_pd(_8,k4));
				  k5 = f(_mm256_add_pd(x0,h2_3),
				         _mm256_fmadd_pd(h81,_mm256_sub_pd(t2,t3),y0));
				  t0 = _mm256_fmadd_pd(_6,k1,
				                   _mm256_fmadd_pd(_36,k2,
						               _mm256_fmadd_pd(_10,k3,
							                   _mm256_mul_pd(_8,k4))));
				  k6 = f(_mm256_add_pd(x0,h4_5),
				         _mm256_fmadd_pd(h75,t0,y0));
				  t1 = _mm256_fmadd_pd(_23,k1,
				                   _mm256_fmsub_pd(_125,k3,_81));
				  t2 = _mm256_fmadd_pd(_125,k6,k5);
				  y0 = _mm256_fmadd_pd(h192,_mm256_mul_pd(t1,t2),y0);
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
			nystrom_step_ymm8r4(__m256(*f)(__m256,
			                                __m256),
					    __m256 y0,
					    __m256 x0,
					    const __m256 h,
					    const int32_t n) {

			    const __m256 _1_25  = _mm256_set1_ps(0.04f);
			    const __m256 h25    = _mm256_mul_ps(_1_25,h);
			    const __m256 _1_4   = _mm256_set1_ps(0.25f);
			    const __m256 h4     = _mm256_mul_ps(_1_4,h);
			    const __m256 _1_81  = _mm256_set1_ps(0.012345679012345679012345679012f);
			    const __m256 h81    = _mm256_mul_ps(_1_81,h);
			    const __m256 _1_75  = _mm256_set1_ps(0.013333333333333333333333333333f);
			    const __m256 h75    = _mm256_mul_ps(_1_45,h);
			    const __m256 _1_3   = _mm256_set1_ps(0.333333333333333333333333333333f);
			    const __m256 h3     = _mm256_mul_ps(_1_3,h);
			    const __m256 _2_5   = _mm256_set1_ps(0.4f);
			    const __m256 h2_5   = _mm256_mul_ps(_2_5,h);
			    const __m256 _2_3   = _mm256_set1_ps(0.666666666666666666666666666667f);
			    const __m256 h2_3   = _mm256_mul_ps(_2_3,h);
			    const __m256 _4_5   = _mm256_set1_ps(0.8f);
			    const __m256 h4_5   = _mm256_set1_ps(_4_5,h);
			    const __m256 _1_192 = _mm256_set1_ps(0.005208333333333333333333333333f);
			    const __m256 h192   = _mm256_mul_ps(_1_192,h);
			    const __m256 _4     = _mm256_set1_ps(4.0f);
			    const __m256 _6     = _mm256_set1_ps(6.0f);
			    const __m256 _12    = _mm256_set1_ps(12.0f);
			    const __m256 _15    = _mm256_set1_ps(15.0f);
			    const __m256 _90    = _mm256_set1_ps(90.0f);
			    const __m256 _50    = _mm256_set1_ps(50.0f);
			    const __m256 _6     = _mm256_set1_ps(6.0f);
			    const __m256 _8     = _mm256_set1_ps(8.0f);
			    const __m256 _36    = _mm256_set1_ps(36.0f);
			    const __m256 _10    = _mm256_set1_ps(10.0f);
			    const __m256 _23    = _mm256_set1_ps(23.0f);
			    const __m256 _125   = _mm256_set1_ps(125.0f);
			    const __m256 _81    = _mm256_set1_ps(81.0f);
			    __m256 k1,k2,k3,k4,k5,k6;
			    __m256 t0,t1,t2,t3,t4;

			    while(--n >= 0) {
                                  k1 = f(x0,y0);
				  k2 = f(_mm256_add_ps(x0,h3),
				         _mm256_fmadd_ps(h3,k1,y0));
				  t0 = _mm256_fmadd_ps(h25,
				                   _mm256_fmadd_ps(_4,k1,
						               _mm256_mul_ps(_6,k2)),y0);
				  k3 = f(_mm256_add_ps(x0,h2_5),t0);
				  t1 = _mm256_fmadd_ps(h4,
				                   _mm256_fmadd_ps(_15,k3,
						               _mm256_sub_ps(k1,
							                 _mm256_mul_ps(_12,k2))),y0);
				  k4 = f(_mm256_add_ps(x0,h),t1);
				  t2 = _mm256_fmadd_ps(_6,k1,
				                   _mm256_mul_ps(_90,k2));
				  t3 = _mm256_fmadd_ps(_50,k3,
				                   _mm256_mul_ps(_8,k4));
				  k5 = f(_mm256_add_ps(x0,h2_3),
				         _mm256_fmadd_ps(h81,_mm256_sub_pd(t2,t3),y0));
				  t0 = _mm256_fmadd_ps(_6,k1,
				                   _mm256_fmadd_ps(_36,k2,
						               _mm256_fmadd_ps(_10,k3,
							                   _mm256_mul_ps(_8,k4))));
				  k6 = f(_mm256_add_ps(x0,h4_5),
				         _mm256_fmadd_ps(h75,t0,y0));
				  t1 = _mm256_fmadd_ps(_23,k1,
				                   _mm256_fmsub_ps(_125,k3,_81));
				  t2 = _mm256_fmadd_ps(_125,k6,k5);
				  y0 = _mm256_fmadd_ps(h192,_mm256_mul_ps(t1,t2),y0);
				  x0 = _mm256_add_ps(x0,h);
			    }
			    return (y0);
                      }


////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Nystrom_Richardson( double (*f)(double, double),       //
//                      double y0, double x0, double h, int number_of_steps,  //
//                                                   int richarson_columns);  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Nystrom method described above       //
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
			nystrom_richardson_ymm4r8(__m256d(*f)(__m256d,
			                                     __m256d),
						 __m256d y0,
						 __m256d x0,
						 const __m256d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m256d richardson[] = {_mm256_set1_pd(1.0/31.0),
			                                                    _mm256_set1_pd(1.0/63.0),
								            _mm256_set1_pd(1.0/127.0),
								            _mm256_set1_pd(1.0/255.0),
								            _mm256_set1_pd(1.0/511.0),
									    _mm256_set1_pd(1.0/1023.0)};
			   constexpr int32_t MAX_COLS = 7;
			   __ATTR_ALIGN__(32) __m256d dt[MAX_COLS];
			   const __m256d _1_2 = _mm512_set1_pd(0.5);
			   __m256d integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = nystrom_step_ymm4r8(f,y0,x0,h_used,n_subints);
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
			nystrom_richardson_ymm8r4(__m256(*f)(__m256,
			                                     __m256),
						 __m256 y0,
						 __m256 x0,
						 const __m256 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(32) const __m256 richardson[] = { _mm256_set1_ps(1.0f/31.0f),
			                                                    _mm256_set1_ps(1.0f/63.0f),
								            _mm256_set1_ps(1.0f/127.0f),
								            _mm256_set1_ps(1.0f/255.0f),
									    _mm256_set1_pd(1.0f/511.0f),
								            _mm256_set1_ps(1.0f/1023.0f)};
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
                                     integral = nystrom_step_ymm8r4(f,y0,x0,h_used,n_subints);
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






#endif /*__GMS_NYSTROM_STEP_AVX_HPP__*/
