

#ifndef __GMS_DERIVATIVE_H__
#define __GMS_DERIVATIVE_H__

namespace file_version {

    const unsigned int gGMS_DERIVATIVE_MAJOR = 1U;
    const unsigned int gGMS_DERIVATIVE_MINOR = 0U;
    const unsigned int gGMS_DERIVATIVE_MICRO = 0U;
    const unsigned int gGMS_DERIVATIVE_FULLVER =
      1000U*gGMS_DERIVATIVE_MAJOR+
      100U*gGMS_DERIVATIVE_MINOR+
      10U*gGMS_DERIVATIVE_MICRO;
    const char * const pgGMS_DERIVATIVE_CREATION_DATE = "10-10-2021 09:44 AM +00200 (SUN 10 OCT 2021 GMT+2)";
    const char * const pgGMS_DERIVATIVE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DERIVATIVE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DERIVATIVE_DESCRIPTION   = "Scalar i.e. not-vectorized derivative implementation."

}

#include <cstdint>
#include "GMS_config.h"


namespace  gms {

             namespace  math {

                       /*
                            Central 5-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
	             
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	               	double
			stencil_5P(double (*f)(double),
			                   double x,
				           double h,
			                   double &err_ro,
				     double &err_tr); 
				     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	               double stencil_5P_optim(double (*f)(double),
			                        double x,
						double h,
						double &abser); 


			  /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	               double stencil_4P( double (*f)(double),
			                   double x,
				           double h,
			                   double &err_ro,
				           double &err_tr); 


			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	              double stencil_4P_optim(double (*f)(double),
			                        double x,
						double h,
						double &abser); 


			

      }

}




#endif /*__GMS_DERIVATIVE_H__*/
