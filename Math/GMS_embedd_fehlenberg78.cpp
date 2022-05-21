
#include <cmath>
#include <algorithm> //min,max
#include "GMS_embedd_fehlenberg78.h"
#include "GMS_cephes.h"

#define ATTEMPTS 12
#define R8_MIN_SCALE_FACTOR 0.125
#define R8_MAX_SCALE_FACTOR 4.0
#define R4_MIN_SCALE_FACTOR 0.125f
#define R4_MAX_SCALE_FACTOR 4.0f


int32_t
embedd_fehlenberg78(double (*f)(double,double),
                       double * __restrict __ATTR_ALIGN__(16) y,
		       double x,
		       double h,
		       double xmax,
		       double * __restrict h_next,
		       double tolerance) {

     constexpr double err_exponent = 0.142857142857142857142857142857;
     __ATTR_ALIGN__(16) double temp_y[2];
     double err;
     double yy;
     int32_t last_ival = 0;
     int32_t i;
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //
     if(__builtin_expect(xmax<x,0) ||
        __builtin_expect(h<=0.0,0)) { return (-2);}
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.
     *h_next = h;
     y[1] = y[0];
     if(__builtin_expect(xmax==x,0)) {  return (0);}
       // Insure that the step size h is not larger than the length of the //
       // integration interval.                                            //
     if(h>(xmax-x)) { h = xmax - x; last_interval = 1;}

        // Redefine the error tolerance to an error tolerance per unit    //
        // length of the integration interval.                            //
     tolerance /= (xmax - x);

        // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
        // maintain an error less than tolerance * (xmax-x) using an     //
        // initial step size of h and initial value: y = y[0]            //
     temp_y[0] = y[0];
     while(x<xmax) {
     scale = 1.0;
     for (i = 0; i < ATTEMPTS; ++i) {
         err = std::fabs(runge_kutta(f,temp_y,x,h));
         if(err== 0.0) { scale = MAX_SCALE_FACTOR; break; }
         yy = (temp_y[0]==0.0) ? tolerance : std::fabs(temp_y[0]);
         scale = 0.8 * std::pow( tolerance * yy / err,err_exponent);
         scale = std::min(std::max(scale,MIN_SCALE_FACTOR),MAX_SCALE_FACTOR);
         if(err<(tolerance * yy)) break;
         h *= scale;
         if(x + h > xmax) h = xmax - x;
         else if( x + h + 0.5 * h > xmax ) h = 0.5 * h;
      }
      if ( i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
      temp_y[0] = temp_y[1];         
      x += h;
      h *= scale;
      *h_next = h;
      if ( last_interval ) break;
      if (  x + h > xmax ) { last_interval = 1; h = xmax - x; }
      else if ( x + h + 0.5 * h > xmax ) h = 0.5 * h;
   }
    y[1] = temp_y[1];
    return 0;
}

__ATTR_ALIGN__(32)
__ATTR_HOT__
static
double
runge_kutta(double (*f)(double,double),
            double * __restrict __ATTR_ALIGN__(16) y,
	    const double x0,
	    const double h) {

    constexpr double c_1_11 = 41.0 / 840.0;
    constexpr double c6 = 34.0 / 105.0;
    constexpr double c_7_8= 9.0 / 35.0;
    constexpr double c_9_10 = 9.0 / 280.0;

    constexpr double a2 = 2.0 / 27.0;
    constexpr double a3 = 1.0 / 9.0;
    constexpr double a4 = 1.0 / 6.0;
    constexpr double a5 = 5.0 / 12.0;
    constexpr double a6 = 1.0 / 2.0;
    constexpr double a7 = 5.0 / 6.0;
    constexpr double a8 = 1.0 / 6.0;
    constexpr double a9 = 2.0 / 3.0;
    constexpr double a10 = 1.0 / 3.0;

    constexpr double b31 = 1.0 / 36.0;
    constexpr double b32 = 3.0 / 36.0;
    constexpr double b41 = 1.0 / 24.0;
    constexpr double b43 = 3.0 / 24.0;
    constexpr double b51 = 20.0 / 48.0;
    constexpr double b53 = -75.0 / 48.0;
    constexpr double b54 = 75.0 / 48.0;
    constexpr double b61 = 1.0 / 20.0;
    constexpr double b64 = 5.0 / 20.0;
    constexpr double b65 = 4.0 / 20.0;
    constexpr double b71 = -25.0 / 108.0;
    constexpr double b74 =  125.0 / 108.0;
    constexpr double b75 = -260.0 / 108.0;
    constexpr double b76 =  250.0 / 108.0;
    constexpr double b81 = 31.0/300.0;
    constexpr double b85 = 61.0/225.0;
    constexpr double b86 = -2.0/9.0;
    constexpr double b87 = 13.0/900.0;
    constexpr double b91 = 2.0;
    constexpr double b94 = -53.0/6.0;
    constexpr double b95 = 704.0 / 45.0;
    constexpr double b96 = -107.0 / 9.0;
    constexpr double b97 = 67.0 / 90.0;
    constexpr double b98 = 3.0;
    constexpr double b10_1 = -91.0 / 108.0;
    constexpr double b10_4 = 23.0 / 108.0;
    constexpr double b10_5 = -976.0 / 135.0;
    constexpr double b10_6 = 311.0 / 54.0;
    constexpr double b10_7 = -19.0 / 60.0;
    constexpr double b10_8 = 17.0 / 6.0;
    constexpr double b10_9 = -1.0 / 12.0;
    constexpr double b11_1 = 2383.0 / 4100.0;
    constexpr double b11_4 = -341.0 / 164.0;
    constexpr double b11_5 = 4496.0 / 1025.0;
    constexpr double b11_6 = -301.0 / 82.0;
    constexpr double b11_7 = 2133.0 / 4100.0;
    constexpr double b11_8 = 45.0 / 82.0;
    constexpr double b11_9 = 45.0 / 164.0;
    constexpr double b11_10 = 18.0 / 41.0;
    constexpr double b12_1 = 3.0 / 205.0;
    constexpr double b12_6 = - 6.0 / 41.0;
    constexpr double b12_7 = - 3.0 / 205.0;
    constexpr double b12_8 = - 3.0 / 41.0;
    constexpr double b12_9 = 3.0 / 41.0;
    constexpr double b12_10 = 6.0 / 41.0;
    constexpr double b13_1 = -1777.0 / 4100.0;
    constexpr double b13_4 = -341.0 / 164.0;
    constexpr double b13_5 = 4496.0 / 1025.0;
    constexpr double b13_6 = -289.0 / 82.0;
    constexpr double b13_7 = 2193.0 / 4100.0;
    constexpr double b13_8 = 51.0 / 82.0;
    constexpr double b13_9 = 33.0 / 164.0;
    constexpr double b13_10 = 12.0 / 41.0;
   
    constexpr double err_factor  = -41.0 / 840.0;
    double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
    double h2_7 = a2 * h;

    k1 = (*f)(x0, *y);
    k2 = (*f)(x0+h2_7, *y + h2_7 * k1);
    k3 = (*f)(x0+a3*h, *y + h * ( b31*k1 + b32*k2) );
    k4 = (*f)(x0+a4*h, *y + h * ( b41*k1 + b43*k3) );
    k5 = (*f)(x0+a5*h, *y + h * ( b51*k1 + b53*k3 + b54*k4) );
    k6 = (*f)(x0+a6*h, *y + h * ( b61*k1 + b64*k4 + b65*k5) );
    k7 = (*f)(x0+a7*h, *y + h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6) );
    k8 = (*f)(x0+a8*h, *y + h * ( b81*k1 + b85*k5 + b86*k6 + b87*k7) );
    k9 = (*f)(x0+a9*h, *y + h * ( b91*k1 + b94*k4 + b95*k5 + b96*k6
                                                          + b97*k7 + b98*k8) );
    k10 = (*f)(x0+a10*h, *y + h * ( b10_1*k1 + b10_4*k4 + b10_5*k5 + b10_6*k6
                                          + b10_7*k7 + b10_8*k8 + b10_9*k9 ) );
    k11 = (*f)(x0+h, *y + h * ( b11_1*k1 + b11_4*k4 + b11_5*k5 + b11_6*k6
                           + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10 ) );
    k12 = (*f)(x0, *y + h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8
                                                 + b12_9*k9 + b12_10 * k10 ) );
    k13 = (*f)(x0+h, *y + h * ( b13_1*k1 + b13_4*k4 + b13_5*k5 + b13_6*k6
                + b13_7*k7 + b13_8*k8 + b13_9*k9 + b13_10*k10 + k12 ) );
    *(y+1) = *y +  h * (c_1_11 * (k1 + k11)  + c6 * k6 + c_7_8 * (k7 + k8) 
                                           + c_9_10 * (k9 + k10) );
    return (err_factor * (k1 + k11 - k12 - k13));
}



int32_t
embedd_fehlenberg78(float (*f)(float,float),
                       float * __restrict __ATTR_ALIGN__(8) y,
		       float x,
		       float h,
		       float xmax,
		       float * __restrict h_next,
		       float tolerance) {

     constexpr float err_exponent = 0.142857142857142857142857142857f;
     __ATTR_ALIGN__(8) float temp_y[2];
     float err;
     float yy;
     int32_t last_ival = 0;
     int32_t i;
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //
     if(__builtin_expect(xmax<x,0) ||
        __builtin_expect(h<=0.0f,0)) { return (-2);}
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.
     *h_next = h;
     y[1] = y[0];
     if(__builtin_expect(xmax==x,0)) {  return (0);}
       // Insure that the step size h is not larger than the length of the //
       // integration interval.                                            //
     if(h>(xmax-x)) { h = xmax - x; last_interval = 1;}

        // Redefine the error tolerance to an error tolerance per unit    //
        // length of the integration interval.                            //
     tolerance /= (xmax - x);

        // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
        // maintain an error less than tolerance * (xmax-x) using an     //
        // initial step size of h and initial value: y = y[0]            //
     temp_y[0] = y[0];
     while(x<xmax) {
     scale = 1.0f;
     for (i = 0; i < ATTEMPTS; ++i) {
         err = std::fabsf(runge_kutta(f,temp_y,x,h));
         if(err== 0.0f) { scale = MAX_SCALE_FACTOR; break; }
         yy = (temp_y[0]==0.0f) ? tolerance : std::fabsf(temp_y[0]);
         scale = 0.8f * cephes_powf( tolerance * yy / err,err_exponent);
         scale = std::min(std::max(scale,MIN_SCALE_FACTOR),MAX_SCALE_FACTOR);
         if(err<(tolerance * yy)) break;
         h *= scale;
         if(x + h > xmax) h = xmax - x;
         else if( x + h + 0.5f * h > xmax ) h = 0.5f * h;
      }
      if( i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
      temp_y[0] = temp_y[1];         
      x += h;
      h *= scale;
      *h_next = h;
      if ( last_interval ) break;
      if (  x + h > xmax ) { last_interval = 1; h = xmax - x; }
      else if ( x + h + 0.5f * h > xmax ) h = 0.5f * h;
   }
    y[1] = temp_y[1];
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Fehlberg's embedded 7th and 8th order methods to     //
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
float
runge_kutta(float (*f)(float,float),
            float * __restrict __ATTR_ALIGN__(16) y,
	    const float x0,
	    const float h) {

    constexpr float c_1_11 = 41.0f / 840.0f;
    constexpr float c6 = 34.0f / 105.0f;
    constexpr float c_7_8= 9.0f / 35.0f;
    constexpr float c_9_10 = 9.0f / 280.0f;

    constexpr float a2 = 2.0f / 27.0f;
    constexpr float a3 = 1.0f / 9.0f;
    constexpr float a4 = 1.0f / 6.0f;
    constexpr float a5 = 5.0f / 12.0f;
    constexpr float a6 = 1.0f / 2.0f;
    constexpr float a7 = 5.0f / 6.0f;
    constexpr float a8 = 1.0f / 6.0f;
    constexpr float a9 = 2.0f / 3.0f;
    constexpr float a10 = 1.0f / 3.0f;

    constexpr float b31 = 1.0f / 36.0f;
    constexpr float b32 = 3.0f / 36.0f;
    constexpr float b41 = 1.0f / 24.0f;
    constexpr float b43 = 3.0f / 24.0f;
    constexpr float b51 = 20.0f / 48.0f;
    constexpr float b53 = -75.0f / 48.0f;
    constexpr float b54 = 75.0f / 48.0f;
    constexpr float b61 = 1.0f / 20.0f;
    constexpr float b64 = 5.0f / 20.0f;
    constexpr float b65 = 4.0f / 20.0f;
    constexpr float b71 = -25.0f / 108.0f;
    constexpr float b74 =  125.0f / 108.0f;
    constexpr float b75 = -260.0f / 108.0f;
    constexpr float b76 =  250.0f / 108.0f;
    constexpr float b81 = 31.0f/300.0f;
    constexpr float b85 = 61.0f/225.0f;
    constexpr float b86 = -2.0f/9.0f;
    constexpr float b87 = 13.0f/900.0f;
    constexpr float b91 = 2.0f;
    constexpr float b94 = -53.0f/6.0f;
    constexpr float b95 = 704.0f / 45.0f;
    constexpr float b96 = -107.0f / 9.0f;
    constexpr float b97 = 67.0f / 90.0f;
    constexpr float b98 = 3.0f;
    constexpr float b10_1 = -91.0f / 108.0f;
    constexpr float b10_4 = 23.0f / 108.0f;
    constexpr float b10_5 = -976.0f / 135.0f;
    constexpr float b10_6 = 311.0f / 54.0f;
    constexpr float b10_7 = -19.0f / 60.0f;
    constexpr float b10_8 = 17.0f / 6.0f;
    constexpr float b10_9 = -1.0f / 12.0f;
    constexpr float b11_1 = 2383.0f / 4100.0f;
    constexpr float b11_4 = -341.0f / 164.0f;
    constexpr float b11_5 = 4496.0f / 1025.0f;
    constexpr float b11_6 = -301.0f / 82.0f;
    constexpr float b11_7 = 2133.0f / 4100.0f;
    constexpr float b11_8 = 45.0f / 82.0f;
    constexpr float b11_9 = 45.0f / 164.0f;
    constexpr float b11_10 = 18.0f / 41.0f;
    constexpr float b12_1 = 3.0f / 205.0f;
    constexpr float b12_6 = - 6.0f / 41.0f;
    constexpr float b12_7 = - 3.0f / 205.0f;
    constexpr float b12_8 = - 3.0f / 41.0f;
    constexpr float b12_9 = 3.0f / 41.0f;
    constexpr float b12_10 = 6.0f / 41.0f;
    constexpr float b13_1 = -1777.0f / 4100.0f;
    constexpr float b13_4 = -341.0 / 164.0;
    constexpr float b13_5 = 4496.0f / 1025.0f;
    constexpr float b13_6 = -289.0f / 82.0f;
    constexpr float b13_7 = 2193.0f / 4100.0f;
    constexpr float b13_8 = 51.0f / 82.0f;
    constexpr float b13_9 = 33.0f / 164.0f;
    constexpr float b13_10 = 12.0f / 41.0f;
   
    constexpr float err_factor  = -41.0f / 840.0f;
    float k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
    float h2_7 = a2 * h;

    k1 = (*f)(x0, *y);
    k2 = (*f)(x0+h2_7, *y + h2_7 * k1);
    k3 = (*f)(x0+a3*h, *y + h * ( b31*k1 + b32*k2) );
    k4 = (*f)(x0+a4*h, *y + h * ( b41*k1 + b43*k3) );
    k5 = (*f)(x0+a5*h, *y + h * ( b51*k1 + b53*k3 + b54*k4) );
    k6 = (*f)(x0+a6*h, *y + h * ( b61*k1 + b64*k4 + b65*k5) );
    k7 = (*f)(x0+a7*h, *y + h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6) );
    k8 = (*f)(x0+a8*h, *y + h * ( b81*k1 + b85*k5 + b86*k6 + b87*k7) );
    k9 = (*f)(x0+a9*h, *y + h * ( b91*k1 + b94*k4 + b95*k5 + b96*k6
                                                          + b97*k7 + b98*k8) );
    k10 = (*f)(x0+a10*h, *y + h * ( b10_1*k1 + b10_4*k4 + b10_5*k5 + b10_6*k6
                                          + b10_7*k7 + b10_8*k8 + b10_9*k9 ) );
    k11 = (*f)(x0+h, *y + h * ( b11_1*k1 + b11_4*k4 + b11_5*k5 + b11_6*k6
                           + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10 ) );
    k12 = (*f)(x0, *y + h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8
                                                 + b12_9*k9 + b12_10 * k10 ) );
    k13 = (*f)(x0+h, *y + h * ( b13_1*k1 + b13_4*k4 + b13_5*k5 + b13_6*k6
                + b13_7*k7 + b13_8*k8 + b13_9*k9 + b13_10*k10 + k12 ) );
    *(y+1) = *y +  h * (c_1_11 * (k1 + k11)  + c6 * k6 + c_7_8 * (k7 + k8) 
                                           + c_9_10 * (k9 + k10) );
    return (err_factor * (k1 + k11 - k12 - k13));
}

