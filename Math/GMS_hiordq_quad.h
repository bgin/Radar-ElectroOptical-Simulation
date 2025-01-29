

#ifndef __GMS_HIORDQ_QUAD_H__
#define __GMS_HIORDQ_QUAD_H__ 081020221251


namespace file_info {


        const unsigned int GMS_HIORDQ_QUAD_MAJOR = 1;

	const unsigned int GMS_HIORDQ_QUAD_MINOR = 1;

	const unsigned int GMS_HIORDQ_QUAD_MICRO = 0;

	const unsigned int GMS_HIORDQ_QUAD_FULLVER = 
		1000U*GMS_HIORDQ_QUAD_MAJOR+100U*GMS_HIORDQ_QUAD_MINOR+10U*GMS_HIORDQ_QUAD_MICRO;

	const char * const GMS_HIORDQ_QUAD_CREATE_DATE = "08-10-2022 12:51 +00200 (SAT 08 OCT 2022 GMT+2)";

	const char * const GMS_HIORDQ_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_HIORDQ_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}

/*
!*****************************************************************************80
!
!! HIORDQ approximates the integral of a function using equally spaced data.
!
!  Discussion:
!
!    The method applies the trapezoidal rule to various subsets of the
!    data, and then applies Richardson extrapolation.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    Alan Kaylor Cline,
!    Department of Computer Science,
!    University of Texas at Austin.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, number of data points.
!
!    Input, real ( kind = 8 ) DELT, the spacing between the X values of the
!    data.  The actual X values are not needed!
!
!    Input, real ( kind = 8 ) Y(NTAB), the Y values of the data.
!
!    Work array, real ( kind = 8 ) WORK(2*(NTAB-1)).  The actual minimum amount
!    of workspace required is two times the number of integer
!    divisors of NTAB-1.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.

*/

#include <cstdint>
#include "GMS_config.h"


namespace gms {

          namespace math {
                
                       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void hiordq(const int32_t ntab,
                                    const double delt,
                                    const double * __restrict __ATTR_ALIGN__(64) y,
                                    const double * __restrict __ATTR_ALIGN__(64) work,
                                    double &result); 


                        
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void hiordq(const int32_t ntab,
                                    const float delt,
                                    const float * __restrict __ATTR_ALIGN__(64) y,
                                    const float * __restrict __ATTR_ALIGN__(64) work,
                                    float &result); 


   
        } // math

} //gms





#endif /*__GMS_HIORDQ_QUAD_H__*/
