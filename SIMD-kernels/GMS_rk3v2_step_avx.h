

#ifndef __GMS_RK3V2_STEP_AVX_H__
#define __GMS_RK3V2_STEP_AVX_H__ 020620221301



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


                    
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               
			__m256d
			rk3v2_step_ymm4r8(__m256d(*f)(__m256d,
			                             __m256d),
					  __m256d y0,
					  __m256d x0,
					  const __m256d h,
					  const int32_t n); 


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               
			__m512
			rk3v2_step_ymm8r4(__m512(*f)(__m512,
			                             __m512),
					  __m512 y0,
					  __m512 x0,
					  const __m512 h,
					  const int32_t n); 

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



                       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	               
			__m256d
			rk3v2_richardson_ymm4r8(__m256d(*f)(__m256d,
			                                     __m256d),
						 __m256d y0,
						 __m256d x0,
						 const __m256d h,
						 const int32_t n,
						 int32_t n_cols);


		        
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	              
			__m512
			rk3v2_richardson_ymm8r4(__m512(*f)(__m512,
			                                     __m512),
						 __m512 y0,
						 __m512 x0,
						 const __m512 h,
						 const int32_t n,
						 int32_t n_cols); 

		      

     }

}







#endif /*__GMS_RK3V2_STEP_AVX_H__*/
