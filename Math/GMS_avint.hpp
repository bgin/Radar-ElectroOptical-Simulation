

#ifndef __GMS_AVINT_HPP__
#define __GMS_AVINT_HPP__ 040920221106


namespace file_info {


        const unsigned int GMS_AVINT_MAJOR = 1;

	const unsigned int GMS_AVINT_MINOR = 1;

	const unsigned int GMS_AVINT_MICRO = 0;

	const unsigned int GMS_AVINT_FULLVER = 
		1000U*GMS_AVINT_MAJOR+100U*GMS_AVINT_MINOR+10U*GMS_AVINT_MICRO;

	const char * const GMS_AVINT_CREATE_DATE = "04-09-2022 11:06 +00200 (SUN 04 SEP 2022 GMT+2)";

	const char * const GMS_AVINT_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_AVINT_AUTHOR = "Data Integrator ported from Fortran version";

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
#include <limits> //Nan value
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
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        float avint(const float * __restrict x,
                                    const float * __restrict y,
                                    const int32_t n,
                                    const float xlo,
                                    const float xup,
				    float &ans,
                                    int32_t &ierr)        {

                          float a , b , c , ca , cb , cc , fl , fr , r3 , &
                                rp5 , slope , sum , syl , syl2 , syl3 , syu , &
                                syu2 , syu3 , term1 , term2 , term3 , x1 , &
                                x12 , x13 , x2 , x23 , x3;
                          int32_t i , inlft , inrt , istart , istop;
                          ierr = 1;
                          ans = 0.0f;
                          // Error check and trivial cases
                          if(xlo==xup) return;
                          if(xlo>xup) {
                             ierr = 2;
                             return;
                          }
                          if(n<2) {
                             ierr = 5;
                             return;
                          }
                          for(i = 1; i < n; ++i) {
                               if(x[i]<=x[i-1]) {
                                  ierr = 4;
                                  return;
                               }
                               if(x[i]>xup) break;
                         }
                         
                         if(n<3) {
                             
                            //! special n=2 case
                            slope = (y[1]-y[0])/(x[1]-x[0]);
                            fl    =  y[0] + slope*(xlo-x[0]);
                            fr    =  y[1] + slope*(xup-x[1]);
                            ans   =  0.5*(fl+fr)*(xup-xlo);
                       
                         }
                         else if(x[n-2]<xlo) {
                            ierr = 3;
                                              
                         }
                         else if(x[2]<=xup) {

                              i = 1;
                              while(true) {
                                
                                 if(x[i]>=xlo) {
                                    inlft = i;
                                    i     = n;
                                    while(true) {
                                        
                                       if(x[i]<=xup) {
                                          inrt = i;
                                          if((inrt-inlft)>=2) {
                                              istart = inlft;
                                              if(inlft==1) istart = 2;
                                              istop = inrt;
                                              if(inrt==n) istop = n-1;
                                              r3  = 3.0f;
                                              rp5 = 0.5f;
                                              sum = 0.0f;
                                              syl = xlo;
                                              syl2= syl*syl;
                                              syl3= syl2*syl;
                                              for(i = istart; i < istop; ++i) {
                                                  x1    = x[i-1];
                                                  x2    = x[i];
                                                  x3    = x[i+1];
                                                  x12   = x1-x2;
                                                  x13   = x1-x3;
                                                  x23   = x2-x3;
                                                  term1 = y[i-1]/(x12*x13);
                                                  term2 = -y[i]/(x12*x23);
                                                  term3 = y[i+1]/(x13*x23);
                                                  a     = term1 + term2 + term3;
                                                  b     = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3;
                                                  c     = x2*x3*term1 + x1*x3*term2 + x1*x2*term3;
                                                  if(i>istart) {
                                                        ca = 0.5f*(a+ca);
                                                        cb = 0.5f*(b+cb);
                                                        cc = 0.5f*(c+cc);
                                                  }
                                                  else {
                                                        ca = a;
                                                        cb = b;
                                                        cc = c;
                                                  }
                                                  syu  = x2;
                                                  syu2 = syu*syu;
                                                  syu3 = syu2*syu;
                                                  sum  = sum + ca*(syu3-syl3)/r3 + cb*rp5*(syu2-syl2) + cc*(syu-syl);
                                                  ca   = a;
                                                  cb   = b;
                                                  cc   = c;
                                                  syl  = syu;
                                                  syl2 = syu2;
                                                  syl3 = syu3;
                                              }//end for
                                               syu = xup;
                                               ans = sum + ca*(syu*syu*syu-syl3)/r3 + cb*rp5*(syu*syu-syl2) + cc*(syu-syl);
                                          }
                                          else {
                                               ierr = 3;
                                          }
                                          return;
                                       } //end if
                                       i -= 1;
                                    } //end while
                                 } // end if
                                  i += 1;
                             }// end while
                         }
                         else {
                              ierr = 3;
                        } //end if
                          
                   }


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        double avint(const double * __restrict x,
                                     const double * __restrict y,
                                     const int32_t n,
                                     const double xlo,
                                     const double xup,
				     double &ans,
                                     int32_t &ierr)        {

                          double a , b , c , ca , cb , cc , fl , fr , r3 , &
                                 rp5 , slope , sum , syl , syl2 , syl3 , syu , &
                                 syu2 , syu3 , term1 , term2 , term3 , x1 , &
                                 x12 , x13 , x2 , x23 , x3;
                          int32_t i , inlft , inrt , istart , istop;
                          ierr = 1;
                          ans = 0.0;
                          // Error check and trivial cases
                          if(xlo==xup) return;
                          if(xlo>xup) {
                             ierr = 2;
                             return;
                          }
                          if(n<2) {
                             ierr = 5;
                             return;
                          }
                          for(i = 1; i < n; ++i) {
                               if(x[i]<=x[i-1]) {
                                  ierr = 4;
                                  return;
                               }
                               if(x[i]>xup) break;
                         }
                         
                         if(n<3) {
                             
                            //! special n=2 case
                            slope = (y[1]-y[0])/(x[1]-x[0]);
                            fl    =  y[0] + slope*(xlo-x[0]);
                            fr    =  y[1] + slope*(xup-x[1]);
                            ans   =  0.5*(fl+fr)*(xup-xlo);
                       
                         }
                         else if(x[n-2]<xlo) {
                            ierr = 3;
                                              
                         }
                         else if(x[2]<=xup) {

                              i = 1;
                              while(true) {
                                
                                 if(x[i]>=xlo) {
                                    inlft = i;
                                    i     = n;
                                    while(true) {
                                        
                                       if(x[i]<=xup) {
                                          inrt = i;
                                          if((inrt-inlft)>=2) {
                                              istart = inlft;
                                              if(inlft==1) istart = 2;
                                              istop = inrt;
                                              if(inrt==n) istop = n-1;
                                              r3  = 3.0;
                                              rp5 = 0.5;
                                              sum = 0.0;
                                              syl = xlo;
                                              syl2= syl*syl;
                                              syl3= syl2*syl;
                                              for(i = istart; i < istop; ++i) {
                                                  x1    = x[i-1];
                                                  x2    = x[i];
                                                  x3    = x[i+1];
                                                  x12   = x1-x2;
                                                  x13   = x1-x3;
                                                  x23   = x2-x3;
                                                  term1 = y[i-1]/(x12*x13);
                                                  term2 = -y[i]/(x12*x23);
                                                  term3 = y[i+1]/(x13*x23);
                                                  a     = term1 + term2 + term3;
                                                  b     = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3;
                                                  c     = x2*x3*term1 + x1*x3*term2 + x1*x2*term3;
                                                  if(i>istart) {
                                                        ca = 0.5*(a+ca);
                                                        cb = 0.5*(b+cb);
                                                        cc = 0.5*(c+cc);
                                                  }
                                                  else {
                                                        ca = a;
                                                        cb = b;
                                                        cc = c;
                                                  }
                                                  syu  = x2;
                                                  syu2 = syu*syu;
                                                  syu3 = syu2*syu;
                                                  sum  = sum + ca*(syu3-syl3)/r3 + cb*rp5*(syu2-syl2) + cc*(syu-syl);
                                                  ca   = a;
                                                  cb   = b;
                                                  cc   = c;
                                                  syl  = syu;
                                                  syl2 = syu2;
                                                  syl3 = syu3;
                                              }//end for
                                               syu = xup;
                                               ans = sum + ca*(syu*syu*syu-syl3)/r3 + cb*rp5*(syu*syu-syl2) + cc*(syu-syl);
                                          }
                                          else {
                                               ierr = 3;
                                          }
                                          return;
                                       } //end if
                                       i -= 1;
                                    } //end while
                                 } // end if
                                  i += 1;
                             }// end while
                         }
                         else {
                              ierr = 3;
                        } //end if
                          
                   }


                   


    }

}
















#endif /*__GMS_AVINT_HPP__*/
