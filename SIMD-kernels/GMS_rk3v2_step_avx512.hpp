

#ifndef __GMS_RK3V2_STEP_AVX512_HPP__
#define __GMS_RK3V2_STEP_AVX512_HPP__ 020620221301



/*
     Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_v2_3.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_RK3V2_STEP_AVX512_MAJOR = 1U;
    const unsigned int GMS_RK3V2_STEP_AVX512_MINOR = 0U;
    const unsigned int GMS_RK3V2_STEP_AVX512_MICRO = 0U;
    const unsigned int GMS_RK3V2_STEP_AVX512_FULLVER =
      1000U*GMS_RK3V2_STEP_AVX512_MAJOR+
      100U*GMS_RK3V2_STEP_AVX512_MINOR+
      10U*GMS_RK3V2_STEP_AVX512_MICRO;
    const char * const GMS_RK3V2_STEP_AVX512_CREATION_DATE = "02-06-2022 13:01 PM +00200 (THR 02 JUN 2022 GMT+2)";
    const char * const GMS_RK3V2_STEPAVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RK3V2_STEP_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RK3V2_STEP_AVX512_DESCRIPTION   = "Vectorized (AVX512) Runge-Kutta order 3 (version 2) step."

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
			__m512d
			rk3v2_step_zmm8r8(__m512d(*f)(__m512d,
			                             __m512d),
					  __m512d y0,
					  __m512d x0,
					  const __m512d h,
					  const int32_t n) {

                           const __m512d _1_3 = _mm512_set1_pd(0.3333333333333333333333333333333);
			   const __m512d _2_3 = _mm512_set1_pd(0.6666666666666666666667);
			   const __m512d _1_4 = _mm512_set1_pd(0.25);
			   const __m512d _3   = _mm512_set1_pd(3.0);
			   const __m512d h13  = _mm512_mul_pd(_1_3,h);
			   const __m512d h23  = _mm512_mul_pd(_2_3,h);
			   const __m512d h4   = _mm512_mul_pd(_1_4,h);
			   __m512d k1,k2;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm512_add_pd(x0,h13),
				        _mm512_fmadd_pd(h13,k1,y0));
				 y0 = _mm512_add_pd(y0,_mm512_mul_pd(h4,
				                                _mm512_fmadd_pd(_3,f(_mm512_add_pd(x0,h23),
								                     _mm512_fmadd_pd(h23,k2,y0)),k1)));
				 x0 = _mm512_add_pd(x0,h);
			   }
			   return (y0);
		      }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			rk3v2_step_zmm16r4(__m512(*f)(__m512,
			                             __m512),
					  __m512 y0,
					  __m512 x0,
					  const __m512 h,
					  const int32_t n) {

                           const __m512 _1_3 = _mm512_set1_ps(0.3333333333333333333333333333333f);
			   const __m512 _2_3 = _mm512_set1_ps(0.6666666666666666666667f);
			   const __m512 _1_4 = _mm512_set1_ps(0.25f);
			   const __m512 _3   = _mm512_set1_ps(3.0f);
			   const __m512 h13  = _mm512_mul_ps(_1_3,h);
			   const __m512 h23  = _mm512_mul_ps(_2_3,h);
			   const __m512 h4   = _mm512_mul_ps(_1_4,h);
			   __m512 k1,k2;

			   while(--n >= 0) {
                                 k1 = f(x0,y0);
				 k2 = f(_mm512_add_ps(x0,h13),
				        _mm512_fmadd_ps(h13,k1,y0));
				 y0 = _mm512_add_ps(y0,_mm512_mul_ps(h4,
				                                _mm512_fmadd_ps(_3,f(_mm512_add_ps(x0,h23),
								                     _mm512_fmadd_ps(h23,k2,y0)),k1)));
				 x0 = _mm512_add_ps(x0,h);
			   }
			   return (y0);
		      }


		      

     }

}







#endif /*__GMS_RK3V2_STEP_AVX512_HPP__*/
