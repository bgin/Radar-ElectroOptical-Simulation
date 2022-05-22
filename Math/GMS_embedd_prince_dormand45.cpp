

#include <cmath>
#include <algorithm> //min,max
#include "GMS_embedd_prince_dormand45.h"
#include "GMS_cephes.h"

#define ATTEMPTS 12
#define R8_MIN_SCALE_FACTOR 0.125
#define R8_MAX_SCALE_FACTOR 4.0
#define R4_MIN_SCALE_FACTOR 0.125f
#define R4_MAX_SCALE_FACTOR 4.0f



int32_t
embedd_prince_dormand45( double (*f)(double,double),
                         double * __restrict __ATTR_ALIGN__(16) y,
		         double x,
		         double h,
		         double xmax, 
		         double * __restrict h_next,
		         double tolerance) {

     __ATTR_ALIGN__(16) double temp_y[2];
     double scale;
     double err;
     double yy;
     int32_t i;
     int32_t last_interval = 0;
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //

     if(__builtin_expect(xmax<x,0) ||
        __builtin_expect(h<= 0.0,0)) return (-2);
   
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.                            //

     *h_next = h;
     y[1] = y[0];
     if(__builtin_expec(xmax==x,0)) return (0); 

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
      while( x < xmax ) {
           scale = 1.0;
           for(i = 0; i < ATTEMPTS; i++) {
                err = std::fabs(runge_kutta(f, temp_y, x, h));
                if(err == 0.0) { scale = R8_MAX_SCALE_FACTOR; break; }
                yy = (temp_y[0] == 0.0) ? tolerance : std::fabs(temp_y[0]);
                scale = 0.8 * std::sqrt(std::sqrt(tolerance * yy/err));
                scale = std::min(std::max(scale,R8_MIN_SCALE_FACTOR),R8_MAX_SCALE_FACTOR);
                if(err < ( tolerance * yy)) break;
                h *= scale;
                if(x + h > xmax) h = xmax - x;
                else if ( x + h + 0.5 * h > xmax ) h = 0.5 * h;
          }
          if(i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
          temp_y[0] = temp_y[1];         
          x += h;
          h *= scale;
          *h_next = h;
          if(last_interval ) break;
          if( x + h > xmax ) { last_interval = 1; h = xmax - x; }
          else if ( x + h + 0.5 * h > xmax ) h = 0.5 * h;
        }
        y[1] = temp_y[1];
        return 0;
}

////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Prince-Dormand's embedded 4th and 5th order methods  //
//     to approximate the solution of the differential equation y'=f(x,y)     //
//     with the initial condition y = y[0] at x = x0.  The value at x + h is  //
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

   constexpr double two_thirds = 2.0 / 3.0;
   constexpr double one_seventwoninths = 1.0 / 729.0;
   constexpr double one_twoninesevenzero = 1.0 / 2970.0;
   constexpr double one_twofivetwozero = 1.0 / 2520.0;
   constexpr double one_fiveninefourzero = 1.0 / 5940.0;

   double k1, k2, k3, k4, k5, k6;
   double h5 = 0.2 * h;

   k1 = (*f)(x0, *y);
   k2 = (*f)(x0+h5, *y + h5 * k1);
   k3 = (*f)(x0+0.3*h, *y + h * ( 0.075 * k1 + 0.225 * k2) );
   k4 = (*f)(x0+0.6*h, *y + h * ( 0.3 * k1 - 0.9 * k2 + 1.2 * k3) );
   k5 = (*f)(x0+two_thirds * h,  *y + one_seventwoninths * h * ( 226.0 * k1 
                                - 675.0 * k2 + 880.0 * k3 + 55.0 * k4) );
   k6 = (*f)(x0+h, *y + one_twoninesevenzero * h * ( - 1991 * k1 + 7425.0 * k2
                               - 2660.0 * k3 - 10010.0 * k4 + 10206.0 * k5) );
   *(y+1) = *y +  one_fiveninefourzero * h * ( 341.0 * k1 + 3800.0 * k3
                                   - 7975.0 * k4 + 9477.0 * k5 + 297.0 * k6 );
   return one_twofivetwozero * (77.0 * k1 - 400.0 * k3 + 1925.0 * k4
                                               - 1701.0 * k5 + 99.0 * k6);
}


int32_t
embedd_prince_dormand45( float (*f)(float,float),
                         float * __restrict __ATTR_ALIGN__(8) y,
		         float x,
		         float h,
		         float xmax, 
		         float * __restrict h_next,
		         float tolerance) {

     __ATTR_ALIGN__(8) float temp_y[2];
     float scale;
     float err;
     float yy;
     int32_t i;
     int32_t last_interval = 0;
      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial enpoint.               //

     if(__builtin_expect(xmax<x,0) ||
        __builtin_expect(h<= 0.0f,0)) return (-2);
   
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.                            //

     *h_next = h;
     y[1] = y[0];
     if(__builtin_expec(xmax==x,0)) return (0); 

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
      while( x < xmax ) {
           scale = 1.0f;
           for(i = 0; i < ATTEMPTS; i++) {
                err = std::fabsf(runge_kutta(f, temp_y, x, h));
                if(err == 0.0) { scale = R4_MAX_SCALE_FACTOR; break; }
                yy = (temp_y[0] == 0.0) ? tolerance : std::fabsf(temp_y[0]);
                scale = 0.8 * cephes_sqrtf(cephes_sqrtf(tolerance * yy/err));
                scale = std::min(std::max(scale,R4_MIN_SCALE_FACTOR),R4_MAX_SCALE_FACTOR);
                if(err < (tolerance * yy)) break;
                h *= scale;
                if(x + h > xmax) h = xmax - x;
                else if(x + h + 0.5f * h > xmax ) h = 0.5f * h;
          }
          if(i >= ATTEMPTS ) { *h_next = h * scale; return -1; };
          temp_y[0] = temp_y[1];         
          x += h;
          h *= scale;
          *h_next = h;
          if(last_interval ) break;
          if(x + h > xmax ) { last_interval = 1; h = xmax - x; }
          else if(x + h + 0.5f * h > xmax ) h = 0.5f * h;
        }
        y[1] = temp_y[1];
        return 0;
}


__ATTR_ALIGN__(32)
__ATTR_HOT__
static
float
runge_kutta(float (*f)(float,float),
            float * __restrict __ATTR_ALIGN__(8) y,
	    const float x0,
	    const float h) {

   constexpr float two_thirds = 2.0f / 3.0f;
   constexpr float one_seventwoninths = 1.0f / 729.0;
   constexpr float one_twoninesevenzero = 1.0f / 2970.0;
   constexpr float one_twofivetwozero = 1.0f / 2520.0f;
   constexpr float one_fiveninefourzero = 1.0f / 5940.0f;

   float k1, k2, k3, k4, k5, k6;
   float h5 = 0.2f * h;

   k1 = (*f)(x0, *y);
   k2 = (*f)(x0+h5, *y + h5 * k1);
   k3 = (*f)(x0+0.3f*h, *y + h * ( 0.075f * k1 + 0.225f * k2) );
   k4 = (*f)(x0+0.6f*h, *y + h * ( 0.3f * k1 - 0.9f * k2 + 1.2f * k3) );
   k5 = (*f)(x0+two_thirds * h,  *y + one_seventwoninths * h * ( 226.0f * k1 
                                - 675.0f * k2 + 880.0f * k3 + 55.0f  * k4) );
   k6 = (*f)(x0+h, *y + one_twoninesevenzero * h * ( -1991.0f * k1 + 7425.0f * k2
                               - 2660.0f * k3 - 10010.0f * k4 + 10206.0f * k5) );
   *(y+1) = *y +  one_fiveninefourzero * h * ( 341.0f * k1 + 3800.0f * k3
                                   - 7975.0f * k4 + 9477.0f * k5 + 297.0f * k6 );
   return one_twofivetwozero * (77.0f * k1 - 400.0f * k3 + 1925.0f * k4
                                               - 1701.0f * k5 + 99.0f * k6);
}
