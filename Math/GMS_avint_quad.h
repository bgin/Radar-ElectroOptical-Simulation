

#ifndef __GMS_AVINT_QUAD_HPP__
#define __GMS_AVINT_QUAD_HPP__ 040920221106


namespace file_info {


        const unsigned int GMS_AVINT_QUAD_MAJOR = 1;

	const unsigned int GMS_AVINT_QUAD_MINOR = 1;

	const unsigned int GMS_AVINT_QUAD_MICRO = 0;

	const unsigned int GMS_AVINT_QUAD_FULLVER = 
		1000U*GMS_AVINT_QUAD_MAJOR+100U*GMS_AVINT_QUAD_MINOR+10U*GMS_AVINT_QUAD_MICRO;

	const char * const GMS_AVINT_QUAD_CREATE_DATE = "04-09-2022 11:06 +00200 (SUN 04 SEP 2022 GMT+2)";

	const char * const GMS_AVINT_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_AVINT_QUAD_AUTHOR = "Data Integrator ported from Fortran version";

}


/*
 !********************************************************************************
!>
!  Integrate a function tabulated at arbitrarily spaced
!  abscissas using overlapping parabolas.
!
!  DAVINT integrates a function tabulated at arbitrarily spaced
!  abscissas.  The limits of integration need not coincide
!  with the tabulated abscissas.
!
!  A method of overlapping parabolas fitted to the data is used
!  provided that there are at least 3 abscissas between the
!  limits of integration.  DAVINT also handles two special cases.
!  If the limits of integration are equal, DAVINT returns a
!  result of zero regardless of the number of tabulated values.
!  If there are only two function values, DAVINT uses the
!  trapezoid rule.
!
!  DAVINT is documented completely in SC-M-69-335
!  Original program from *Numerical Integration* by Davis & Rabinowitz
!  Adaptation and modifications by Rondall E Jones.
!
!### References
!  * R. E. Jones, Approximate integrator of functions
!    tabulated at arbitrarily spaced abscissas,
!    Report SC-M-69-335, Sandia Laboratories, 1969.
!
!### Author
!  * Jones, R. E., (SNLA)
!
!### REVISION HISTORY
!  * 690901  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jacob Williams, Jan 2022 : modernized this procedure. added quad-precision coefficients.
*/

#include <cstdint>
#include "GMS_config.h"


namespace gms {


      namespace math {

/*
     real(kind=dp),dimension(:),intent(in) :: x !! array of abscissas, which must be in increasing order.
    real(kind=dp),dimension(:),intent(in) :: y !! array of function values. i.e., `y(i)=func(x(i))`
    integer(kind=i4),intent(in) :: n !! The integer(kind=i4) number of function values supplied.
                            !! `N >= 2` unless `XLO = XUP`.
    real(kind=dp),intent(in) :: xlo !! lower limit of integration
    real(kind=dp),intent(in) :: xup !! upper limit of integration.  Must have `XLO <= XUP`
    real(kind=dp),intent(out) :: ans !! computed approximate value of integral
    integer(kind=i4),intent(out) :: ierr !! A status code:
                                !!
                                !! * Normal Code
                                !!    * =1 Means the requested integration was performed.
                                !! * Abnormal Codes
                                !!    * =2 Means `XUP` was less than `XLO`.
                                !!    * =3 Means the number of `X(I)` between `XLO` and `XUP`
                                !!      (inclusive) was less than 3 and neither of the two
                                !!      special cases described in the abstract occurred.
                                !!      No integration was performed.
                                !!    * =4 Means the restriction `X(I+1)>X(I)` was violated.
                                !!    * =5 Means the number `N` of function values was < 2.
                                !!
                                !! ANS is set to zero if `IERR` = 2, 3, 4, or 5.
*/
                   
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        float avint(const float * __restrict x,
                                    const float * __restrict y,
                                    const int32_t n,
                                    const float xlo,
                                    const float xup,
                                    int32_t &ierr);        
                                    
                                    

                     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        double avint(const double * __restrict x,
                                     const double * __restrict y,
                                     const int32_t n,
                                     const double xlo,
                                     const double xup,
                                     int32_t &ierr);       

                   


    }

}
















#endif /*__GMS_AVINT_QUAD_H__*/
