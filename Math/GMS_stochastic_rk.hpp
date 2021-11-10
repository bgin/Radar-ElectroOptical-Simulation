
#ifndef __GMS_STOCHASTIC_RK_HPP__
#define __GMS_STOCHASTIC_RK_HPP__ 101120211408



namespace file_info {

      const unsigned int gGMS_STOCHASTIC_RK_PD_MAJOR = 1;
      const unsigned int gGMS_STOCHASTIC_RK_PD_MINOR = 0;
      const unsigned int gGMS_STOCHASTIC_RK_PD_MICRO = 0;
      const unsigned int gGMS_STOCHASTIC_RK_PD_FULLVER =
        1000*gGMS_STOCHASTIC_RK_PD_MAJOR+100*gGMS_STOCHASTIC_RK_PD_MINOR+
	10*gGMS_STOCHASTIC_RK_PD_MICRO;
      const char * const pgGMS_STOCHASTIC_RK_PD_BUILD_DATE = __DATE__":"__TIME__;
      const char * const pgGMS_STOCHASTIC_RK_PD_CREATION_DATE = "10-11-2021 14:08  +00200 (WED 10 NOV 2021 GMT+2)";
      const char * const pgGMS_STOCHASTIC_RK__PD_DESCRIPTION   = "Stochastic Runge-Kutte AVX512 scalar."


}



/*
     Modified:
    07 July 2010
  Author:
    John Burkardt
  Modified:  
    
     Bernard Gingold on 05-11-2021 14:16  +00200 (FRI 05 NOV 2021 GMT+2)
     Original implementation slightly modified.
     Removed J. Burkardt pseudorandom (scalar) generators.
     Added gcc function attributes.
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
*/



#include <cstdint>
#include <cmath>
#include "GMS_config.h"


 /*
         The Runge-Kutta scheme is first-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
     Parameters:

    Input, double X, the value at the current time.

    Input, double T, the current time.

    Input, double H, the time step.

    Input, double Q, the spectral density of the input white noise.

    Input, double FI ( double X ), the name of the deterministic
    right hand side function.

    Input, double GI ( double X ), the name of the stochastic
    right hand side function.

    Input double rand, random value passed from the callee routine.

    Output, double RK1_TI_STEP, the value at time T+H.
*/

                                    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            static inline
                                    double rk1_ti_step_scalar(const double x,
				                              const double t,
							      const double h,
							      const double q,
							      const double rand,
							      double (*fi) (const double),
							      double (*gi) (const double)) {

                                            constexpr double a21 = 0.0;
					    constexpr double q1  = 0.0;
                                            double k1            = 0.0;
                                            double w1            = 0.0;
                                            double step          = 0.0;
                                            w1 = rand  * std::sqrt(q1 * q / h);
                                            k1 = h * fi( x1) + h * gi( x1) * w1;
                                            step = x1 + a21 * k1;
					    return (step);
				   }

				   

                                    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            static inline
                                    double rk2_ti_step_scalar( const double x,
				                              const double t,
							      const double h,
							      const double q,
							      const double rand1,
							      const double rand2,
							      double (*fi) (const double),
							      double (*gi) (const double)) {

                                            constexpr double a21 = 1.0;
                                            constexpr double a31 = 0.5;
                                            constexpr double a32 = 0.5;
					    constexpr double q1  = 2.0;
                                            constexpr double q2  = q1;
					    const     double qh  = q/h;
                                            double k1 = 0.0;
                                            double k2 = 0.0;
                                            double w1 = 0.0;
                                            double w2 = 0.0;
                                            double x2 = 0.0;
                                            double step = 0.0;
                                            w1 = rand1 * std::sqrt(q1 * qh);
                                            k1 = h * fi( x) + h * gi( x) * w1;
                                            x2 = x1 + a21 * k1;
                                            w2 = rand2 * std::sqrt(q2 * qh);
                                            k2 = h * fi( x2) + h * gi( x2) * w2;
                                            step = x1 + a31 * k1 + a32 * k2;
					    return (step);
				  }



				    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            static inline
                                    double rk3_ti_step_scalar( const double x,
				                              const double t,
							      const double h,
							      const double q,
							      const double rand1,
							      const double rand2,
							      const double rand3,
							      double (*fi) (const double),
							      double (*gi) (const double)) {

                                         constexpr double a21 =  1.52880952525675;
                                         constexpr double a31 =  0.0;
                                         constexpr double a32 =  0.51578733443615;
                                         constexpr double a41 =  0.53289582961739;
                                         constexpr double a42 =  0.25574324768195;
                                         constexpr double a43 =  0.21136092270067;
					 const double     qh  =  q/h;
                                         double k1 = 0.0;
                                         double k2 = 0.0;
                                         double k3 = 0.0;
                                         constexpr double q1  =  1.87653936176981;
                                         constexpr double q2  =  3.91017166264989;
                                         constexpr double q3  =  4.73124353935667;
                                         double w1 = 0.0;
                                         double w2 = 0.0;
                                         double w3 = 0.0;
                                         double x2 = 0.0;
                                         double x3 = 0.0;
                                         double step = 0.0;
                                         w1 = rand1 * std::sqrt(q1 * qh);
                                         k1 = h * fi(x) + h * gi(x) * w1;
                                         x2 = x1 + a21 * k1;
                                         w2 = rand2 * std::sqrt(q2 * qh);
                                         k2 = h * fi( x2) + h * gi( x2) * w2;
                                         x3 = x1 + a31 * k1 + a32 * k2;
                                         w3 = rand3 * std::sqrt(q3 * qh);
                                         k3 = h * fi( x3) + h * gi( x3) * w3;
                                         step = x1 + a41 * k1 + a42 * k2 + a43 * k3;
                                         return (step);
				 }


				 
                                    __ATTR_ALWAYS_INLINE__
				    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            static inline
                                    double rk4_ti_step_scalar(const double x,
				                              const double t,
							      const double h,
							      const double q,
							      const double rand1,
							      const double rand2,
							      const double rand3,
							      const double rand4,
							      double (*fi) (const double),
							      double (*gi) (const double)) {

					       constexpr  double a21 =  2.71644396264860;
                                               constexpr  double a31 = -6.95653259006152;
                                               constexpr  double a32 =  0.78313689457981;
                                               constexpr  double a41 =  0.0;
                                               constexpr  double a42 =  0.48257353309214;
                                               constexpr  double a43 =  0.26171080165848;
                                               constexpr  double a51 =  0.47012396888046;
                                               constexpr  double a52 =  0.36597075368373;
                                               constexpr  double a53 =  0.08906615686702;
                                               constexpr  double a54 =  0.07483912056879;
					       const double qh       = q/h;
                                               double k1 = 0.0;
                                               double k2 = 0.0;
                                               double k3 = 0.0;
                                               double k4 = 0.0;
                                               constexpr double q1 = 2.12709852335625;
                                               constexpr double q2 = 2.73245878238737;
                                               constexpr double q3 = 11.22760917474960;
                                               constexpr double q4 = 13.36199560336697;
                                               double w1;
                                               double w2;
                                               double w3;
                                               double w4;
                                               double x2;
                                               double x3;
                                               double x4;
                                               double step;
                                               w1 = rand1 * std::sqrt(q1 * qh);
                                               k1 = h * fi(x) + h * gi(x) * w1;
                                               x2 = x1 + a21 * k1;
                                               w2 = rand2 * std::sqrt(q2 * qh);
                                               k2 = h * fi(x2) + h * gi(x2) * w2;
                                               x3 = x1 + a31 * k1 + a32 * k2;
                                               w3 = rand3 * std::sqrt(q3 * qh);
                                               k3 = h * fi(x3) + h * gi(x3) * w3;
                                               x4 = x1 + a41 * k1 + a42 * k2;
                                               w4 = rand4 * std::sqrt(q4 * qh);
                                               k4 = h * fi(x4) + h * gi(x4) * w4;
                                               step = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;
					       return (step);

				 }



#endif /*__GMS_STOCHASTIC_RK_HPP__*/
