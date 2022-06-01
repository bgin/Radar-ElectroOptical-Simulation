
#ifndef __GMS_RK4_STEP_AVX512_HPP__
#define __GMS_RK4_STEP_AVX512_HPP__

/*
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_RK4_STEP_AVX512_MAJOR = 1U;
    const unsigned int GMS_RK4_STEP_AVX512_MINOR = 0U;
    const unsigned int GMS_RK4_STEP_AVX512_MICRO = 0U;
    const unsigned int GMS_RK4_STEP_AVX512_FULLVER =
      1000U*GMS_RK4_STEP_AVX512_MAJOR+
      100U*GMS_RK4_STEP_AVX512_MINOR+
      10U*GMS_RK4_STEP_AVX512_MICRO;
    const char * const GMS_RK4_STEP_AVX512_CREATION_DATE = "01-06-2022 08:27 PM +00200 (WED 01 JUN 2022 GMT+2)";
    const char * const GMS_RK4_STEPAVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RK4_STEP_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RK4_STEP_AVX512_DESCRIPTION   = "Vectorized (AVX512) Runge-Kutta order 4 step."

}


#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

        namespace math {

/*
////////////////////////////////////////////////////////////////////////////////
//  double Runge_Kutta( double (*f)(double, double), double y0, double x0,    //
//                                          double h, int number_of_steps );  //
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
//                                                                            //
*/

	                __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d
		        rk4_step_zmm8r8(__m512d(*f)(const __m512d,
			                            const __m512d),
				        __m512d y0,
				        __m512d x0,
					const __m512d h,
					const int32_t n) {

                           const __m512d _1_6 = _mm512_set1_pd(0.1666666666666666666667);
			   const __m512d _0_5 = _mm512_set1_pd(0.5);
			   const __m512d h2   = _mm512_mul_pd(_0_5,h);
			   __m512d k1,k2,k3,k4,tmp;
			   while(--n >= 0) { // for 'n' attempts do:
                                 k1 = _mm512_mul_pd(h,f(x0,y0));
				 k2 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,h2),
				                                  _mm512_fmadd_pd(_0_5,k1,y0)));
				 tmp= _mm512_add_pd(k1,k2);
				 k3 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,h2),
				                                  _mm512_fmadd_pd(_0_5,k2,y0)));
				 tmp= _mm512_add_pd(tmp,_mm512_add_pd(k2,k3));
				 x0 = _mm512_add_pd(x0,h);
				 k4 = _mm512_mul_pd(h,f(x0,_mm512_add_pd(y0,k3)));
				 tmp= _mm512_add_pd(tmp,_mm512_add_pd(k3,k4));
				 y0 = _mm512_fmadd_pd(y0,_1_6,tmp);
			  }
			  return (y0);
		     }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
		        rk4_step_zmm16r4(__m512(*f)(const __m512,
			                            const __m512),
				         __m512 y0,
				         __m512 x0,
					 const __m512 h,
					 const int32_t n) {

                           const __m512 _1_6 = _mm512_set1_ps(0.1666666666666666666667f);
			   const __m512 _0_5 = _mm512_set1_ps(0.5f);
			   const __m512 h2   = _mm512_mul_ps(_0_5,h);
			   __m512 k1,k2,k3,k4,tmp;
			   while(--n >= 0) { // for 'n' attempts do:
                                 k1 = _mm512_mul_ps(h,f(x0,y0));
				 k2 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,h2),
				                                  _mm512_fmadd_ps(_0_5,k1,y0)));
				 tmp= _mm512_add_ps(k1,k2);
				 k3 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,h2),
				                                  _mm512_fmadd_ps(_0_5,k2,y0)));
				 tmp= _mm512_add_ps(tmp,_mm512_add_ps(k2,k3));
				 x0 = _mm512_add_ps(x0,h);
				 k4 = _mm512_mul_ps(h,f(x0,_mm512_add_ps(y0,k3)));
				 tmp= _mm512_add_ps(tmp,_mm512_add_ps(k3,k4));
				 y0 = _mm512_fmadd_ps(y0,_1_6,tmp);
			  }
			  return (y0);
		     }
			

   } // math

} //gms









#endif /*__GMS_RK4_STEP_AVX512_HPP__*/
