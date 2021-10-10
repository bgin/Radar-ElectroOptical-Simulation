

#ifndef __GMS_DERIVATIVE_HPP__
#define __GMS_DERIVATIVE_HPP__

namespace file_version {

    const unsigned int gGMS_DERIVATIVE_MAJOR = 1U;
    const unsigned int gGMS_DERIVATIVE_MINOR = 0U;
    const unsigned int gGMS_DERIVATIVE_MICRO = 0U;
    const unsigned int gGMS_DERIVATIVE_FULLVER =
      1000U*gGMS_DERIVATIVE_MAJOR+
      100U*gGMS_DERIVATIVE_MINOR+
      10U*gGMS_DERIVATIVE_MICRO;
    const char * const pgGMS_DERIVATIVE_CREATION_DATE = "10-10-2021 09:44 AM +00200 (SUN 10 OCT 2021 GMT+2)";
    const char * const pgGMS_SETV512_AVX_UNROLL16X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DERIVATIVE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DERIVATIVE_DESCRIPTION   = "Scalar i.e. not-vectorized derivative implementation."

}

#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include "GMS_config.h"


namespace  gms {

             namespace  math {

                       /*
                            Central 5-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
	                __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
	                double
			stencil_5P(double (*f)(double),
			                   double x,
				           double h,
			                   double &err_ro,
				     double &err_tr) {
                            
                             constexpr double n0_5  = 0.5;
			     constexpr double n1_3  = 0.333333333333333333333333;
			     constexpr double n4_3  = 1.3333333333333333333333;
			     constexpr double n2    = 2.0;
			     double p1;
			     double p2;
			     double p1h;
			     double p2h;
			     double t0;
			     double t1;
			     double t2;
			     double t3;
			     double t4;
			     double tmp0;
			     double tmp1;
			     p1  = f(x-h);
			     p2  = f(x+h);
			     t0  = n0_5*(p1-p2);
			     p1h = f(x-h*n0_5);
			     p2h = f(x+h*n0_5);
			     t1  = n4_3*(p2h-p1h)-n3_0*t0;
			     t2  = (std::fabs(p2)+std::fabs(p1))*
			           std::numeric_limits<double>::eps();
			     tmp0 = t1/h;
			     t3  = n2*(std::fabs(p2h)+std::fabs(p1h))*
			           std::numeric_limits<double>::eps()+t2;
			     tmp1 = t2/h;
			     t4  = std::max(tmp0,tmp1)*(std::fabs(x)/h)*
			           std::numeric_limits<double>::eps();
			     err_ro = std::fabs(t3/h)*t4;
			     err_tr = std::fabs((t1-t0)/h);
			     return (t1/h);
			}


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
	                double stencil_5P_optim(double (*f)(double),
			                        double x,
						double h,
						double &abser) {
                               constexpr double n2 = 2.0;
			       constexpr double n1_3 = 0.333333333333333333333333333;
			       constexpr double n4 = 4.0;
			       double x0   = 0.0;
			       double err_ro = 0.0;
                               double err_tr = 0.0;
			       double err;
			       x0 = stencil_5P(f,x,h,err_ro,err_tr);
			       err = err_ro+err_tr;
			       const bool b1 = err_ro>0.0 && err_tr>0.0;
			       const bool b2 = err_ro<err_tr;
			       if(b1 && b2) {
                                  double x02 = 0.0;
				  double err_ro2 = 0.0;
				  double err_tr2 = 0.0;
				  double err2;
				  double h2;
				  h2 = h*std::pow(err_ro/(err_tr+err_tr),n1_3);
				  x02=stencil_5P(f,x,h,err_ro2,err_tr2);
				  err2 = err_ro2+err_tr2;
				  const bool b3 = err2<err;
				  const bool b4 = std::fabs(x02-x0)<n4*err;
				  if(b3 && b4) {
                                     x0  = x02;
				     err = err2;
				  }
			       }
			       abser = err;
			       return (x0);
			}


			  /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
	                double stencil_4P( double (*f)(double),
			                   double x,
				           double h,
			                   double &err_ro,
				           double &err_tr) {

                             constexpr double n1_4  = 0.25;
			     constexpr double n1_2  = 0.5;
			     constexpr double n3_4  = 0.75;
			     constexpr double n22_3 = 7.3333333333333333333333;
			     constexpr double n62_3 = 20.6666666666666666666667;
			     constexpr double n52_3 = 17.3333333333333333333333;
			     constexpr double n2    = 2.0;
			     double dydx = 0.0;
			     double p1;
			     double p2;
			     double p3;
			     double p4;
			     double t0;
			     double t1;
			     double t2;
			     double t3;
			     p1 = f(x+h*n1_4);
			     p2 = f(x+h*n1_2);
			     p3 = f(x+n3_4*h);
			     p4 = f(x+h);
			     t0 = n2*(p4-p2);
			     t1 = n22_3*(p4-p3)-n62_3*(p3-p2)+
			          n52_3*(p2-p1);
			     t2 = 41.34*(std::fabs(p4)+std::fabs(p3)+
			                 std::fabs(p2)+std::fabs(p1))*
					 std::numeric_limits<double>::eps();
			     t3 = std::max(std::fabs(t0/h),std::fabs(t1/h))*
			          std::fabs(x/h)*std::numeric_limits<double>::eps();
			     dydx = t1/h;
			     err_tr = std::fabs(t2-t0/h);
			     err_ro = std::fabs(t2/h)+t3;
			     return (dydx);
			}


			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
	                double stencil_4P_optim(double (*f)(double),
			                        double x,
						double h,
						double &abser) {
			       constexpr double n1_2 = 0.5;
			       constexpr double n4   = 4.0;
                               double x0   = 0.0;
			       double err_ro = 0.0;
                               double err_tr = 0.0;
			       double err;
			       x0 = stencil_5P(f,x,h,err_ro,err_tr);
			       err = err_ro+err_tr;
			       const bool b1 = err_ro>0.0 && err_tr>0.0;
			       const bool b2 = err_ro<err_tr;
			       if(b1 && b2) {
                                  double x02 = 0.0;
				  double err_ro2 = 0.0;
				  double err_tr2 = 0.0;
				  double err2;
				  double h2;
				  h2 = h * std::pow(err_ro/err_tr,n1_2);
				  x02=stencil_4P(f,x,h,err_ro2,err_tr2);
				  err2 = err_ro2+err_tr2;
				  const bool b3 = err2<err;
				  const bool b4 = std::fabs(x02-x0)<n4*err;
				  if(b3 && b4) {
                                     x0  = x02;
				     err = err2;
				  }
			       }
			       abser = err;
			       return (x0);
			}


			

      }

}




#endif /*__GMS_DERIVATIVE_HPP__*/
