
#ifndef __GMS_PLINT_QUAD_H__
#define __GMS_PLINT_QUAD_H__


namespace file_info {


        const unsigned int GMS_PLINT_QUAD_MAJOR = 1;

	const unsigned int GMS_PLINT_QUAD_MINOR = 1;

	const unsigned int GMS_PLINT_QUAD_MICRO = 0;

	const unsigned int GMS_PLINT_QUAD_FULLVER = 
		1000U*GMS_PLINT_QUAD_MAJOR+100U*GMS_PLINT_QUAD_MINOR+10U*GMS_PLINT_QUAD_MICRO;

	const char * const GMS_PLINT_QUAD_CREATE_DATE = "09-10-2022 16:05 +00200 (SUN 09 OCT 2022 GMT+2)";

	const char * const GMS_PLINT_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_PLINT_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}


/*
   !*****************************************************************************80
!
!! PLINT approximates the integral of unequally spaced data.
!
!  Discussion:
!
!    The method uses piecewise linear interpolation.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 2.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), the function values, 
!    FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.       
*/

#include <cstdint>
#include <limits>
#include "GMS_config.h"


namespace  gms {

          namespace  math {
                  
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        void plint(const int32_t ntab,
                                   double * __restrict __ATTR_ALIGN__(64) xtab,
                                   double * __restrict __ATTR_ALIGN__(64) ftab,
                                   const double a,
                                   const double b,
                                   double &result); 
                                   

                     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			 void plint(const int32_t ntab,
                                   float * __restrict __ATTR_ALIGN__(64) xtab,
                                   float * __restrict __ATTR_ALIGN__(64) ftab,
                                   const float a,
                                   const float b,
                                   double &result); 















 
      } // math

} // gms














#endif /*__GMS_PLINT_QUAD_H__*/
