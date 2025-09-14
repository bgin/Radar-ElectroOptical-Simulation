
#include <cmath>
#include <algorithm> //min,max
#include "GMS_embedd_fehlenberg34.h"
#include "GMS_cephes.h"

#define ATTEMPTS 12
#define R8_MIN_SCALE_FACTOR 0.125
#define R8_MAX_SCALE_FACTOR 4.0
#define R8_ORDER            3.0
#define R4_MIN_SCALE_FACTOR 0.125f
#define R4_MAX_SCALE_FACTOR 4.0f
#define R4_ORDER            3.0f

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Fehlberg method is an adaptive procedure for approxi-  //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y) five//
//     times per step using embedded third order and fourth order Runge-Kutta //
//     estimates to estimate the not only the solution but also the error.    //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * ( 229 / 1470 * k1 + 1125 / 1813 * k3           //
//                              + 13718 / 81585 * k4 + 1 / 18 * k5 )          //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+2h/7, y[i] + 2/7h*k1 ),                                   //
//     k3 = f( x[i]+7h/15, y[i]+h*(77/900 k1 + 343/900 k2) ),                 //
//     k4 = f( x[i]+35h/38, y[i]+h*(805/1444 k1 - 77175/54872 k2              //
//                                                     + 97125/54872 k3) ),   //
//     k5 = f( x[i]+h, y[i]+h*(79/490 k1 + 2175/3626 k3 + 2166/9065 k4)),     //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h * ( - 888 k1 + 3375 k3 - 11552 k4 + 9065 k5 ) / 163170      //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/3     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////


int32_t
embedd_fehlenberg34(double (*f)(double,double),
                       double * __restrict __ATTR_ALIGN__(16) y,
		       double x,
		       double h,
		       double xmax,
		       double * __restrict h_next,
		       double tolerance) {

        constexpr double err_exponent = 0.33333333333333333333333333333333;
	double temp_y[2];
        double scale;
        double err;
        double yy;
        int32_t i;
        int32_t last_interval = 0;
  
   
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //

        if(__builtin_expect(xmax<x,0) ||
	   __builtin_expect(h<=0.0,0)) { return (-2);}
   
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.                            //

       *h_next = h;
       y[1] = y[0];
       if(__builtin_expect(xmax==x,0)) {return (0);} 

       // Insure that the step size h is not larger than the length of the //
       // integration interval.                                            //
  
       h = min(h, xmax - x);

        // Redefine the error tolerance to an error tolerance per unit    //
        // length of the integration interval.                            //

       tolerance /= (xmax - x);

        // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
        // maintain an error less than tolerance * (xmax-x) using an     //
        // initial step size of h and initial value: y = y[0]            //

       temp_y[0] = y[0];
       while ( x < xmax ) {
          scale = 1.0;
          for(i = 0; i < ATTEMPTS; ++i) {
                 err = std::fabs(runge_kutta(f, temp_y, x, h));
                 if(err==0.0) {scale = R8_MAX_SCALE_FACTOR; break; }
                 yy = (temp_y[0] == 0.0) ? tolerance : std::fabs(temp_y[0]);
                 scale = 0.8 * std::pow(tolerance * yy / err, err_exponent);
                 scale = std::min(std::max(scale,R8_MIN_SCALE_FACTOR),R8_MAX_SCALE_FACTOR);
                 if( err < ( tolerance * yy ) ) break;
                 h *= scale;
                 if(x + h > xmax) h = xmax - x;
                 else if ( x + h + 0.5 * h > xmax ) h = 0.5 * h;
          }
          if(i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
          temp_y[0] = temp_y[1];         
          x += h;
          h *= scale;
          *h_next = h;
          if(last_interval) break;
          if(x + h > xmax) { last_interval = 1; h = xmax - x; }
          else if(x + h + 0.5 * h > xmax) h = 0.5 * h;
       }
       y[1] = temp_y[1];
       return 0;

}

////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Fehlberg's embedded 3rd and 4th order methods to     //
//     approximate the solution of the differential equation y'=f(x,y) with   //
//     the initial condition y = y[0] at x = x0.  The value at x + h is       //
//     returned in y[1].  The function returns err / h ( the absolute error   //
//     per step size ).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y[0]).                   //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at x + h.                              //
//     double x   Initial value of x.                                         //
//     double h   Step size                                                   //
//                                                                            //
//  Return Values:                                                            //
//     This routine returns the err / h.  The solution of y(x) at x + h is    //
//     returned in y[1].                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

__ATTR_ALIGN__(32)
__ATTR_HOT__
static
double
runge_kutta(double (*f)(double,double),
            double * __restrict __ATTR_ALIGN__(16) y,
	    const double x0,
	    const double h) {

     constexpr double a2 = 2.0 / 7.0;
     constexpr double a3 = 7.0 / 15.0;
     constexpr double a4 = 35.0 / 38.0;

     constexpr double b31 = 77.0 / 900.0;
     constexpr double b32 = 343.0 / 900.0;
     constexpr double b41 = 805.0 / 1444.0;
     constexpr double b42 = -77175.0 / 54872.0;
     constexpr double b43 = 97125.0 / 54872.0;
     constexpr double b51 = 79.0 / 490.0;
     constexpr double b53 = 2175.0 / 3626.0;
     constexpr double b54 = 2166.0 / 9065.0;

     constexpr double c1 = 229.0 / 1470.0;
     constexpr double c3 = 1125.0 / 1813.0;
     constexpr double c4 = 13718.0 / 81585.0;
     constexpr double c5 = 1.0 / 18.0;

     constexpr double d1 = -888.0 / 163170.0;
     constexpr double d3 = 3375.0 / 163170.0;
     constexpr double d4 = -11552.0 / 163170.0;
     constexpr double d5 = 9065.0 / 163170.0;
     double k1, k2, k3, k4, k5;
     double h2 = a2 * h;

     k1 = (*f)(x0, *y);
     k2 = (*f)(x0+h2, *y + h2 * k1);
     k3 = (*f)(x0+a3*h, *y + h * ( b31*k1 + b32*k2) );
     k4 = (*f)(x0+a4*h, *y + h * ( b41*k1 + b42*k2 + b43*k3) );
     k5 = (*f)(x0+h,  *y + h * ( b51*k1 + b53*k3 + b54*k4) );
     *(y+1) = *y +  h * (c1*k1 + c3*k3 + c4*k4 + c5*k5);
     return  (d1*k1 + d3*k3 + d4*k4 + d5*k5);
}


int32_t
embedd_fehlenberg34(float (*f)(float,float),
                      float * __restrict __ATTR_ALIGN__(16) y,
		      float x,
		      float h,
		      float xmax,
		      float * __restrict h_next,
		      float tolerance) {

        constexpr float err_exponent = 0.33333333333333333333333333333333f;
	__ATTR_ALIGN__(8) float temp_y[2];
        float scale;
        float err;
        float yy;
        int32_t i;
        int32_t last_interval = 0;
  
   
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //

        if(__builtin_expect(xmax<x,0) ||
	   __builtin_expect(h<=0.0f,0)) { return (-2);}
   
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.                            //

       *h_next = h;
       y[1] = y[0];
       if(__builtin_expect(xmax==x,0)) {return (0);} 

       // Insure that the step size h is not larger than the length of the //
       // integration interval.                                            //
  
       h = min(h, xmax - x);

        // Redefine the error tolerance to an error tolerance per unit    //
        // length of the integration interval.                            //

       tolerance /= (xmax - x);

        // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
        // maintain an error less than tolerance * (xmax-x) using an     //
        // initial step size of h and initial value: y = y[0]            //

       temp_y[0] = y[0];
       while ( x < xmax ) {
          scale = 1.0;
          for(i = 0; i < ATTEMPTS; ++i) {
                 err = std::fabsf(runge_kutta(f, temp_y, x, h));
                 if(err==0.0f) {scale = R4_MAX_SCALE_FACTOR; break; }
                 yy = (temp_y[0] == 0.0f) ? tolerance : std::fabsf(temp_y[0]);
                 scale = 0.8f * cephes_powf(tolerance * yy / err, err_exponent);
                 scale = std::min(std::max(scale,R4_MIN_SCALE_FACTOR),R4_MAX_SCALE_FACTOR);
                 if( err < ( tolerance * yy ) ) break;
                 h *= scale;
                 if(x + h > xmax) h = xmax - x;
                 else if ( x + h + 0.5f * h > xmax ) h = 0.5f * h;
          }
          if(i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
          temp_y[0] = temp_y[1];         
          x += h;
          h *= scale;
          *h_next = h;
          if(last_interval) break;
          if(x + h > xmax) { last_interval = 1; h = xmax - x; }
          else if(x + h + 0.5f * h > xmax) h = 0.5f * h;
       }
       y[1] = temp_y[1];
       return 0;

}


__ATTR_ALIGN__(32)
__ATTR_HOT__
static
float
runge_kutta(float (*f)(float,float),
            dfloat * __restrict __ATTR_ALIGN__(16) y,
	    const float x0,
	    const float h) {

     constexpr float a2 = 2.0f / 7.0f;
     constexpr float a3 = 7.0f / 15.0f;
     constexpr float a4 = 35.0f / 38.0f;

     constexpr float b31 = 77.0f / 900.0f;
     constexpr float b32 = 343.0f / 900.0f;
     constexpr float b41 = 805.0f / 1444.0f;
     constexpr float b42 = -77175.0f / 54872.0f;
     constexpr float b43 = 97125.0f / 54872.0f;
     constexpr float b51 = 79.0f / 490.0f;
     constexpr float b53 = 2175.0f / 3626.0f;
     constexpr float b54 = 2166.0f / 9065.0f;

     constexpr float c1 = 229.0f / 1470.0f;
     constexpr float c3 = 1125.0f / 1813.0f;
     constexpr float c4 = 13718.0f / 81585.0f;
     constexpr float c5 = 1.0f / 18.0f;

     constexpr float d1 = -888.0f / 163170.0f;
     constexpr float d3 = 3375.0f / 163170.0f;
     constexpr float d4 = -11552.0f / 163170.0f;
     constexpr float d5 = 9065.0f / 163170.0f;
     float k1, k2, k3, k4, k5;
     float h2 = a2 * h;

     k1 = (*f)(x0, *y);
     k2 = (*f)(x0+h2, *y + h2 * k1);
     k3 = (*f)(x0+a3*h, *y + h * ( b31*k1 + b32*k2) );
     k4 = (*f)(x0+a4*h, *y + h * ( b41*k1 + b42*k2 + b43*k3) );
     k5 = (*f)(x0+h,  *y + h * ( b51*k1 + b53*k3 + b54*k4) );
     *(y+1) = *y +  h * (c1*k1 + c3*k3 + c4*k4 + c5*k5);
     return  (d1*k1 + d3*k3 + d4*k4 + d5*k5);
}
