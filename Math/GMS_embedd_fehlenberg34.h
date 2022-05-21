

#ifndef __GMS_EMBEDD_FEHLENBERG34_H__
#define __GMS_EMBEDD_FEHLENBERG34_H__ 210520221351

/*
   Adapted from following implementation: 
   http://www.mymathlib.com/diffeq/embedded_runge_kutta/embedded_fehlberg_3_4.html
   
*/

namespace file_info {


        const unsigned int GMS_EMBEDD_FEHLENBERG34_MAJOR = 1;

	const unsigned int GMS_EMBEDD_FEHLENBERG34_MINOR = 1;

	const unsigned int GMS_EMBEDD_FEHLENBERG34_MICRO = 0;

	const unsigned int GMS_EMBEDD_FEHLENBERG34_FULLVER = 
		1000U*GMS_EMBEDD_FEHLENBERG34_MAJOR+100U*GMS_EMBEDD_FEHLENBERG34_MINOR+10U*GMS_EMBEDD_FEHLENBERG34_MICRO;

	const char * const GMS_EMBEDD_FEHLENBERG34_CREATE_DATE = "21-05-2022 13:51 +00200 (SAT 21 MAY 2022 GMT+2)";

	const char * const GMS_EMBEDD_FEHLENBERG78_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_EMBEDD_FEHLENBERG78_AUTHOR = "Adapted from: http://www.mymathlib.com/diffeq/embedded_runge_kutta";

}



#include <cstdint>
#include "GMS_config.h"

////////////////////////////////////////////////////////////////////////////////
// int Embedded_Fehlberg_3_4( double (*f)(double, double), double y[],        //
//       double x, double h, double xmax, double *h_next, double tolerance )  //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation y'=f(x,y) with the      //
//     initial condition y(x) = y[0].  The value at xmax is returned in y[1]. //
//     The function returns 0 if successful or -1 if it fails.                //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y0) corresponding to the //
//                initial condition y(x0) = y0.                               //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at xmax.                               //
//     double x   The initial value of x.                                     //
//     double h   Initial step size.                                          //
//     double xmax The endpoint of x.                                         //
//     double *h_next   A pointer to the estimated step size for successive   //
//                      calls to Runge_Kutta_Fehlberg.                        //
//     double tolerance The tolerance of y(xmax), i.e. a solution is sought   //
//                so that the relative error < tolerance.                     //
//                                                                            //
//  Return Values:                                                            //
//     0   The solution of y' = f(x,y) from x to xmax is stored y[1] and      //
//         h_next has the value to the next size to try.                      //
//    -1   The solution of y' = f(x,y) from x to xmax failed.                 //
//    -2   Failed because either xmax < x or the step size h <= 0.            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


int32_t
embedd_fehlenberg34(double (*)(double,double),
                       double * __restrict,
		       double,
		       double,
		       double,
		       double * __restrict,
		       double) __ATTR_ALIGN__(32)
		               __ATTR_HOT__;

	 
int32_t
embedd_fehlenberg34(float (*)(float,float),
                       float * __restrict,
		       float,
		       float,
		       float,
		       float * __restrict,
		       float) __ATTR_ALIGN__(32)
		               __ATTR_HOT__;









#endif /*__GMS_EMBEDD_FEHLENBERG34_H__*/
