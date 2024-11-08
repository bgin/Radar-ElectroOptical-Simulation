




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
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_butcher.html
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/





#include "GMS_butcher_step_avx512.h"
#include "GMS_simd_utils.hpp"


////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Butcher( double (*f)(double, double), double y0,       //
//                               double x0, double h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Butcher method described above to    //
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


                    
			__m512d
			gms::math::butcher_step_zmm8r8(__m512d(*f)(__m512d,
			                             __m512d),
					    __m512d y0,
					    __m512d x0,
					    const __m512d h,
					    const int32_t n) {

                            const __m512d _1_3  = _mm512_set1_pd(0.3333333333333333333333333333333333);
			    const __m512d h3    = _mm512_mul_pd(_1_3,h);
			    const __m512d _2_3  = _mm512_set1_pd(0.666666666666666666666666666667);
			    const __m512d h2_3  = _mm512_mul_pd(_2_3,h);
			    const __m512d _1_2  = _mm512_set1_pd(0.5);
			    const __m512d h2    = _mm512_mul_pd(_1_2,h);
			    const __m512d _1_12 = _mm512_set1_pd(0.083333333333333333333333333333);
			    const __m512d h12   = _mm512_mul_pd(_1_12,h);
			    const __m512d _1_16 = _mm512_set1_pd(0.0625);
			    const __m512d h16   = _mm512_mul_pd(_1_16,h);
			    const __m512d _1_8  = _mm512_set1_pd(0.125);
			    const __m512d h8    = _mm512_mul_pd(_1_8,h);
			    const __m512d _1_44 = _mm512_set1_pd(0.022727272727272727272727272727);
			    const __m512d h44   = _mm512_mul_pd(_1_44,h);
			    const __m512d _1_120= _mm512_set1_pd(0.008333333333333333333333333333);
			    const __m512d h120  = _mm512_mul_pd(_1_120,h);
			    const __m512d _4    = _mm512_set1_pd(4.0);
			    const __m512d _18   = _mm512_set1_pd(18.0);
			    const __m512d _3    = _mm512_set1_pd(3.0);
			    const __m512d _6    = _mm512_set1_pd(6.0);
			    const __m512d _9    = _mm512_set1_pd(9.0);
			    const __m512d _36   = _mm512_set1_pd(36.0);
			    const __m512d _63   = _mm512_set1_pd(63.0);
			    const __m512d _72   = _mm512_set1_pd(72.0);
			    const __m512d _64   = _mm512_set1_pd(64.0);
			    const __m512d _11   = _mm512_set1_pd(11.0);
			    const __m512d _81   = _mm512_set1_pd(81.0);
			    const __m512d _32   = _mm512_set1_pd(32.0);
			    __m512d k1,k2,k3,k4,k5,k6,k7;
			    __m512d t0,t1,t2,t3,t4;

			    while(--n >= 0) {
                                  k1 = f(x0,y0);
				  k2 = f(_mm512_add_pd(x0,h3),
				         _mm512_fmadd_pd(h3,k1,y0));
				  k3 = f(_mm512_add_pd(x0,h2_3),
				         _mm512_fmadd_pd(h2_3,k2,y0));
				  k4 = f(_mm512_add_pd(x0,h3),
				         _mm512_fmadd_pd(h12,
					             _mm512_add_pd(k1,
						               _mm512_fmsub_pd(_4,k2,k3)),y0));
				  t0 = _mm512_sub_pd(_mm512_fmadd_pd(_18,k2,zmm8r8_negate(k1)),
				                     _mm512_fmsub_pd(_3,k3,_mm512_mul_pd(_6,k4)));
				  k5 = f(_mm512_add_pd(x0,h2),
				         _mm512_fmadd_pd(h16,t0,y0));
				  t1 = _mm512_fmsub_pd(_9,k2,
				                   _mm512_fmsub_pd(_3,k3,
						               _mm512_fmadd_pd(_6,k4,
							                   _mm512_mul_pd(_4,k5))));
				  k6 = f(_mm512_add_pd(x0,h2),
				         _mm512_fmadd_pd(h8,t1,y0));
				  t2 = _mm512_fmsub_pd(_9,k1,
				                   _mm512_fmadd_pd(_36,k2,
						               _mm512_fmadd_pd(_63,k3,
							                   _mm512_fmsub_pd(_72,k4,
									               _mm512_mul_pd(_64,k6)))));
				  k7 = f(_mm512_add_pd(x0,h),
				         _mm512_fmadd_pd(h44,t2,y0));
				  t3 = _mm512_mul_pd(_11,_mm512_add_pd(k1,k7));
				  t4 = _mm512_mul_pd(_81,_mm512_add_pd(k3,k4));
				  t0 = _mm512_mul_pd(_32,_mm512_add_pd(k5,k6));
				  y0 = _mm512_fmadd_pd(h120,
				                   _mm512_add_pd(t3,
						             _mm512_sub_pd(t4,t0)),y0);
				  x0 = _mm512_add_pd(x0,h);
					             
			    }
			    return (y0);
		     }


		     
			__m512
			gms::math::butcher_step_zmm16r4(__m512(*f)(__m512,
			                             __m512),
					    __m512 y0,
					    __m512 x0,
					    const __m512 h,
					    const int32_t n) {

                            const __m512 _1_3  = _mm512_set1_ps(0.3333333333333333333333333333333333f);
			    const __m512 h3    = _mm512_mul_ps(_1_3,h);
			    const __m512 _2_3  = _mm512_set1_ps(0.666666666666666666666666666667f);
			    const __m512 h2_3  = _mm512_mul_ps(_2_3,h);
			    const __m512 _1_2  = _mm512_set1_ps(0.5f);
			    const __m512 h2    = _mm512_mul_ps(_1_2,h);
			    const __m512 _1_12 = _mm512_set1_ps(0.083333333333333333333333333333f);
			    const __m512 h12   = _mm512_mul_ps(_1_12,h);
			    const __m512 _1_16 = _mm512_set1_ps(0.0625f);
			    const __m512 h16   = _mm512_mul_ps(_1_16,h);
			    const __m512 _1_8  = _mm512_set1_ps(0.125f);
			    const __m512 h8    = _mm512_mul_ps(_1_8,h);
			    const __m512 _1_44 = _mm512_set1_ps(0.022727272727272727272727272727f);
			    const __m512 h44   = _mm512_mul_ps(_1_44,h);
			    const __m512 _1_120= _mm512_set1_ps(0.008333333333333333333333333333f);
			    const __m512 h120  = _mm512_mul_ps(_1_120,h);
			    const __m512 _4    = _mm512_set1_ps(4.0f);
			    const __m512 _18   = _mm512_set1_ps(18.0f);
			    const __m512 _3    = _mm512_set1_ps(3.0f);
			    const __m512 _6    = _mm512_set1_ps(6.0f);
			    const __m512 _9    = _mm512_set1_ps(9.0f);
			    const __m512 _36   = _mm512_set1_ps(36.0f);
			    const __m512 _63   = _mm512_set1_ps(63.0f);
			    const __m512 _72   = _mm512_set1_ps(72.0f);
			    const __m512 _64   = _mm512_set1_ps(64.0f);
			    const __m512 _11   = _mm512_set1_ps(11.0f);
			    const __m512 _81   = _mm512_set1_ps(81.0f);
			    const __m512 _32   = _mm512_set1_ps(32.0f);
			    __m512 k1,k2,k3,k4,k5,k6,k7;
			    __m512 t0,t1,t2,t3,t4;

			    while(--n >= 0) {
                                  k1 = f(x0,y0);
				  k2 = f(_mm512_add_ps(x0,h3),
				         _mm512_fmadd_ps(h3,k1,y0));
				  k3 = f(_mm512_add_ps(x0,h2_3),
				         _mm512_fmadd_ps(h2_3,k2,y0));
				  k4 = f(_mm512_add_ps(x0,h3),
				         _mm512_fmadd_ps(h12,
					             _mm512_add_ps(k1,
						               _mm512_fmsub_ps(_4,k2,k3)),y0));
				  t0 = _mm512_sub_ps(_mm512_fmadd_ps(_18,k2,zmm16r4_negate(k1)),
				                     _mm512_fmsub_ps(_3,k3,_mm512_mul_ps(_6,k4)));
				  k5 = f(_mm512_add_ps(x0,h2),
				         _mm512_fmadd_ps(h16,t0,y0));
				  t1 = _mm512_fmsub_ps(_9,k2,
				                   _mm512_fmsub_ps(_3,k3,
						               _mm512_fmadd_ps(_6,k4,
							                   _mm512_mul_ps(_4,k5))));
				  k6 = f(_mm512_add_ps(x0,h2),
				         _mm512_fmadd_ps(h8,t1,y0));
				  t2 = _mm512_fmsub_ps(_9,k1,
				                   _mm512_fmadd_ps(_36,k2,
						               _mm512_fmadd_ps(_63,k3,
							                   _mm512_fmsub_ps(_72,k4,
									               _mm512_mul_ps(_64,k6)))));
				  k7 = f(_mm512_add_ps(x0,h),
				         _mm512_fmadd_ps(h44,t2,y0));
				  t3 = _mm512_mul_ps(_11,_mm512_add_ps(k1,k7));
				  t4 = _mm512_mul_ps(_81,_mm512_add_ps(k3,k4));
				  t0 = _mm512_mul_ps(_32,_mm512_add_ps(k5,k6));
				  y0 = _mm512_fmadd_ps(h120,
				                   _mm512_add_ps(t3,
						             _mm512_sub_ps(t4,t0)),y0);
				  x0 = _mm512_add_ps(x0,h);
					             
			    }
			    return (y0);
		     }

////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta_Butcher_Richardson( double (*f)(double, double),       //
//                      double y0, double x0, double h, int number_of_steps,  //
//                                                   int richarson_columns);  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta-Butcher method described above       //
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


                  
                   
			__m512d
			gms::math::butcher_richardson_zmm8r8(__m512d(*f)(__m512d,
			                                     __m512d),
						 __m512d y0,
						 __m512d x0,
						 const __m512d h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(64) const __m512d richardson[] = {_mm512_set1_pd(1.0/63.0),
			                                                    _mm512_set1_pd(1.0/127.0),
								            _mm512_set1_pd(1.0/255.0),
								            _mm512_set1_pd(1.0/511.0),
								            _mm512_set1_pd(1.0/1023.0)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(64) __m512d dt[MAX_COLS];
			   const __m512d _1_2 = _mm512_set1_pd(0.5);
			   __m512d integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = butcher_step_zmm8r8(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm512_sub_pd(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm512_fmadd_pd(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm512_mul_pd(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm512_add_pd(x0,h);
			   }
			   return (y0);
		    }


		     
			__m512
			gms::math::butcher_richardson_zmm16r4(__m512(*f)(__m512,
			                                     __m512),
						 __m512 y0,
						 __m512 x0,
						 const __m512 h,
						 const int32_t n,
						 int32_t n_cols) {

			   __ATTR_ALIGN__(64) const __m512 richardson[] = { _mm512_set1_ps(1.0f/63.0f),
			                                                    _mm512_set1_ps(1.0f/127.0f),
								            _mm512_set1_ps(1.0f/255.0f),
								            _mm512_set1_ps(1.0f/511.0f),
								            _mm512_set1_ps(1.0f/1023.0f)};
			   constexpr int32_t MAX_COLS = 6;
			   __ATTR_ALIGN__(64) __m512 dt[MAX_COLS];
			   const __m512 _1_2 = _mm512_set1_ps(0.5f);
			   __m512 integral,delta,h_used;
			   int32_t n_subints;
			   n_cols = std::max(1,std::min(MAX_COLS,n_cols));
			   while(--n >= 0) {
                                 h_used = h;
				 n_subints = 1;
				 #pragma loop_count min(1),max(6),avg(3)
				 for(int32_t j = 0; j < n_cols; ++j) {
                                     integral = butcher_step_zmm16r4(f,y0,x0,h_used,n_subints);
				     for(int32_t k = 0; k < j; ++k) {
                                         delta = _mm512_sub_ps(integral,dt[k]);
					 dt[k] = integral;
					 integral = _mm512_fmadd_ps(richardson[k],delta,integral);
				     }
				     dt[j] = integral;
				     h_used = _mm512_mul_ps(h_used,_1_2);
				     n_subints += n_subints;
				 }
				 y0 = integral;
				 x0 = _mm512_add_ps(x0,h);
			   }
			   return (y0);
		    }



  
