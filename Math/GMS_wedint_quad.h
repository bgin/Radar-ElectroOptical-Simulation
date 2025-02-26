

#ifndef __GMS_WEDINT_QUAD_H__
#define __GMS_WEDINT_QUAD_H__ 061020221649


namespace file_info {


        const unsigned int GMS_WEDINT_QUAD_MAJOR = 1;

	const unsigned int GMS_WEDINT_QUAD_MINOR = 1;

	const unsigned int GMS_WEDINT_QUAD_MICRO = 0;

	const unsigned int GMS_WEDINT_QUAD_FULLVER = 
		1000U*GMS_WEDINT_QUAD_MAJOR+100U*GMS_WEDINT_QUAD_MINOR+10U*GMS_WEDINT_QUAD_MICRO;

	const char * const GMS_WEDINT_QUAD_CREATE_DATE = "06-10-2022 16:49 +00200 (THR 06 OCT 2022 GMT+2)";

	const char * const GMS_WEDINT_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_WEDINT_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}

/*
 *****************************************************************************80
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, is the number of data points.  
!    (NTAB-1) must be divisible by 6.
!
!    Input, real ( kind = 8 ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated data values.
!
!    Output, real ( kind = 8 ) RESULT, is the approximation to the integral.   
*/

#include <cstdint>
#include "GMS_config.h"


namespace  gms {

         namespace  math {
              
                       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void wedint(const int32_t ntab,
                                    const double h,
                                    const double * __restrict __ATTR_ALIGN__(64) ftab,
                                    double &result); 


                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void wedint(const int32_t ntab,
                                    const float h,
                                    const float * __restrict __ATTR_ALIGN__(64) ftab,
                                    double &result); 
     }

}










#endif /*__GMS_WEDINT_QUAD_H__*/
