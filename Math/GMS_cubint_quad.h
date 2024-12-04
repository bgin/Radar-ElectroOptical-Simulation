
#ifndef __GMS_CUBINT_QUAD_H__
#define __GMS_CUBINT_QUAD_H__ 061020220904


namespace file_info {


        const unsigned int GMS_CUBINT_QUAD_MAJOR = 1;

	const unsigned int GMS_CUBINT_QUAD_MINOR = 1;

	const unsigned int GMS_CUBINT_QUAD_MICRO = 0;

	const unsigned int GMS_CUBINT_QUAD_FULLVER = 
		1000U*GMS_CUBINT_QUAD_MAJOR+100U*GMS_CUBINT_QUAD_MINOR+10U*GMS_CUBINT_QUAD_MICRO;

	const char * const GMS_CUBINT_QUAD_CREATE_DATE = "06-10-2022 09:04 +00200 (THR 06 OCT 2022 GMT+2)";

	const char * const GMS_CUBINT_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_CUBINT_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}

/*
   !*****************************************************************************80
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      Integral ( XTAB(IB) <= X <= XTAB(IA) ) F(X) DX
!
!    The routine estimates the error in integration.
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
!    Philip Gill, GF Miller,
!    An algorithm for the integration of unequally spaced data,
!    The Computer Journal, 
!    Number 15, Number 1, 1972, pages 80-83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, integer ( kind = 4 ) IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer ( kind = 4 ) IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the error in
!    integration.
*/


#include <cstdint>
#include "GMS_config.h"

namespace  gms {

        namespace  math {

                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                void cubint(const int32_t ntab,
                                    double * __restrict __ATTR_ALIGN__(64) xtab,
                                    double * __restrict __ATTR_ALIGN__(64) ftab,
                                    const int32_t ia,
                                    const int32_t ib,
                                    double &result,
                                    double &error); 
                                    

                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		       void cubint(const int32_t ntab,
                                    float * __restrict __ATTR_ALIGN__(64) xtab,
                                    float * __restrict __ATTR_ALIGN__(64) ftab,
                                    const int32_t ia,
                                    const int32_t ib,
                                    float &result,
                                    float &error); 
  
                  
       }
}






















#endif /*__GMS_CUBINT_QUAD_H__*/
