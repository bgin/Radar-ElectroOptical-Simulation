

#ifndef __GMS_EMBEDD_PRINCE_DORMAND45_H__
#define __GMS_EMBEDD_PRINCE_DORMAND45_H__


/*
   Adapted from following implementation: 
   http://www.mymathlib.com/diffeq/embedded_runge_kutta
   
*/

namespace file_info {


        const unsigned int GMS_EMBED_PRINCE_DORMAND45_MAJOR = 1;

	const unsigned int GMS_EMBEDD_PRINCE_DORMAND45_MINOR = 1;

	const unsigned int GMS_EMBEDD_PRINCE_DORMAND45_MICRO = 0;

	const unsigned int GMS_EMBEDD_PRINCE_DORMAND45_FULLVER = 
		1000U*GMS_EMBEDD_PRINCE_DORMAND45_MAJOR+100U*GMS_EMBEDD_PRINCE_DORMAND45_MINOR+10U*GMS_EMBEDD_PRINCE_DORMAND45_MICRO;

	const char * const GMS_EMBEDD_PRINCE_DORMAND45_CREATE_DATE = "22-05-2022 12:44 +00200 (SUN 22 MAY 2022 GMT+2)";

	const char * const GMS_EMBEDD__PRINCE_DORMAND45_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_EMBEDD__PRINCE_DORMAND45_AUTHOR = "Adapted from: http://www.mymathlib.com/diffeq/embedded_runge_kutta";

}



#include <cstdint>
#include "GMS_config.h"


////////////////////////////////////////////////////////////////////////////////
// File: embedded_prince_dormand_v1_4_5.c                                     //
// Routines:                                                                  //
//    Embedded_Prince_Dormand_v1_4_5                                          //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     This Runge-Kutta-Prince-Dormand method is an adaptive procedure for    //
//     approximating the solution of the differential equation y'(x) = f(x,y) //
//     with initial condition y(x0) = c.  This implementation evaluates       //
//     f(x,y) six times per step using embedded fourth order and fifth order  //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * ( 31 / 540 * k1 + 190 / 297 * k3               //
//                           - 145 / 108 * k4 + 351 / 220 * k5 + 1/20 * k6 )  //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/5, y[i] + h*k1/5 ),                                     //
//     k3 = f( x[i]+3h/10, y[i]+(h/40)*(3 k1 + 9 k2) ),                       //
//     k4 = f( x[i]+3h/5, y[i]+(h/10)*(3 k1 - 9 k2 + 12 k3) ),                //
//     k5 = f( x[i]+2h/3, y[i]+(h/729)*(226 k1 - 675 k2 + 880 k3 + 55 k4) )   //
//     k6 = f( x[i]+h, y[i]+(h/2970)*(-1991 k1 + 7425 k2 - 2660 k3            //
//                                                  - 10010 k4 + 10206 k5) )  //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( 77 k1 - 400 k3 + 1925 k4 - 1701 k5 + 99 k6 ) / 2520       //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/4     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////


int32_t
embedd_prince_dormand45(double (*)(double,double),
                         double * __restrict,
		         double,
		         double,
		         double,
		         double * __restrict,
		         double) __ATTR_ALIGN__(32)
		               __ATTR_HOT__;

	 
int32_t
embedd_prince_dormand45(float (*)(float,float),
                         float * __restrict,
		         float,
		         float,
		         float,
		         float * __restrict,
		         float) __ATTR_ALIGN__(32)
		               __ATTR_HOT__;








#endif /*__GMS_EMBEDD_PRINCE_DORMAND145_H__*/
