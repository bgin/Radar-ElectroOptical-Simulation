
#ifndef __GMS_PLINT_QUAD_HPP__
#define __GMS_PLINT_QUAD_HPP__


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
                        __ATTR_ALIGN__(32)
			static inline
                        void plint(const int32_t ntab,
                                   double * __restrict __ATTR_ALIGN__(64) xtab,
                                   double * __restrict __ATTR_ALIGN__(64) ftab,
                                   const double a,
                                   const double b,
                                   double &result) {

                          if(__builtin_expect(ntab<2,0) || 
                             __builtin_expect(a==b,0)) {
                             result = std::numeric_limits<double>::quiet_NaN();
                             return;
                          } 
                          for(int32_t i = 1; i != ntab; ++i) {
                              if(xtab[i]<=xtab[i-1]) {
                                 result = std::numeric_limits<double>::quiet_NaN();
                                 return;
                              }
                          }
                          double fa,fb,slope,syl;
                          int32_t i,ihi,ilo,ind;
                          
                          if(b<a) {
                             syl = b;
                             b   = a;
                             a   = syl
                             ind = -1;
                          }
                          else {
                             syl = a;
                             ind = 1;
                          }
                        
                          ilo = ntab+1;
                          for(i = 0; i != ntab; ++i) {
                              if(a<=xtab[i]) {
                                 ilo = i;
                                 break;
                              }
                          }
                          ihi = 0;
                          for(i = ntab; i != 0; --i) {
                               if(xtab[i]<=b) {
                                  ihi = 1;
                                  break;
                               }
                          }

                          if(ihi==0) {
                             slope  = (ftab[1]-ftab[0])/(xtab[1]-xtab[0]);
                             fa     = ftab[0]+slope*(a-xtab[0]);
                             fb     = ftab[0]+slope*(b-xtab[0]);
                             result = 0.5*(b-a)*(fa+fb);
                          }
                          else if(ilo==ntab+1) {
                             slope  = (ftab[ntab]-ftab[ntab-1])/(xtab[ntab]-xtab[ntab-1]);
                             fa     = ftab[ntab-1]+slope*(a-xtab[ntab-1]);
                             fb     = ftab[ntab-1]+slope*(b-xtab[ntab-1]);
                             result = 0.5*(b-a)*(fa+fb);
                          }
                          else if(ihi+1==ilo) {
                             slope  = (ftab[ilo]-ftab[ihi])/(xtab[ilo]-xtab[ihi]);
                             fa     = ftab[ihi]+slope*(a-xtab[ihi]);
                             fb     = ftab[ihi]+slope*(b-xtab[ihi]);
                             result = 0.5*(b-a)*(fa+fb);
                          }
                          else {
                             result = 0.0;

                             for(i = ilo; i != ihi-1; ++i) {
                                 result = result+0.5*(xtab[i+1]-xtab[i])
                                          * (ftab[i]+ftab[i+1]);
                             }

                             if(ilo==1) {
                                 slope  = (ftab[1]-ftab[0])/(xtab[1]-xtab[0]);
                                 fa     = ftab[0]+slope*(a-xtab[0]);
                                 result = result+0.5*(xtab[ilo]-a)*(fa+ftab[ilo]);
                             }
                             else {
                                 slope  = (ftab[ilo]-ftab[ilo-1])/(xtab[ilo]-xtab[ilo-1]);
                                 fa     = ftab[ilo-1]+slope*(a-xtab[ilo-1]);
                                 result = result+0.5*(xtab[ilo]-a)*(fa+ftab[ilo]);
                             }
 
                            if(ihi==ntab) {
                                 slope  = (ftab[ntab]-ftab[ntab-1])/(xtab[ntab]-xtab[ntab-1] ;
                                 fb     = ftab[ntab-1]+slope*(b-xtab[ntab-1]);
                                 result = result+0.5*(b-xtab[ntab])*(fb+ftab[ntab]);
                            }
                            else {
                                 slope  = (ftab[ihi+1]-ftab[ihi])/(xtab[ihi+1]-xtab[ihi]);
                                 fb     = ftab[ihi]+slope*(b-xtab[ihi]);
                                 result = result+0.5*(b-xtab[ihi])*(fb+ftab[ihi]);
                            }
                         }

                         if(ind!=1) {
                             ind = 1;
                             syl = b;
                             b = a;
                             a = syl;
                             result = -result;
                          }
                          
                     } 


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void plint(const int32_t ntab,
                                   float * __restrict __ATTR_ALIGN__(64) xtab,
                                   float * __restrict __ATTR_ALIGN__(64) ftab,
                                   const float a,
                                   const float b,
                                   double &result) {

                          if(__builtin_expect(ntab<2,0) || 
                             __builtin_expect(a==b,0)) {
                             result = std::numeric_limits<float>::quiet_NaN();
                             return;
                          } 
                          for(int32_t i = 1; i != ntab; ++i) {
                              if(xtab[i]<=xtab[i-1]) {
                                 result = std::numeric_limits<float>::quiet_NaN();
                                 return;
                              }
                          }
                          float fa,fb,slope,syl;
                          int32_t i,ihi,ilo,ind;
                          
                          if(b<a) {
                             syl = b;
                             b   = a;
                             a   = syl
                             ind = -1;
                          }
                          else {
                             syl = a;
                             ind = 1;
                          }
                        
                          ilo = ntab+1;
                          for(i = 0; i != ntab; ++i) {
                              if(a<=xtab[i]) {
                                 ilo = i;
                                 break;
                              }
                          }
                          ihi = 0;
                          for(i = ntab; i != 0; --i) {
                               if(xtab[i]<=b) {
                                  ihi = 1;
                                  break;
                               }
                          }

                          if(ihi==0) {
                             slope  = (ftab[1]-ftab[0])/(xtab[1]-xtab[0]);
                             fa     = ftab[0]+slope*(a-xtab[0]);
                             fb     = ftab[0]+slope*(b-xtab[0]);
                             result = 0.5f*(b-a)*(fa+fb);
                          }
                          else if(ilo==ntab+1) {
                             slope  = (ftab[ntab]-ftab[ntab-1])/(xtab[ntab]-xtab[ntab-1]);
                             fa     = ftab[ntab-1]+slope*(a-xtab[ntab-1]);
                             fb     = ftab[ntab-1]+slope*(b-xtab[ntab-1]);
                             result = 0.5f*(b-a)*(fa+fb);
                          }
                          else if(ihi+1==ilo) {
                             slope  = (ftab[ilo]-ftab[ihi])/(xtab[ilo]-xtab[ihi]);
                             fa     = ftab[ihi]+slope*(a-xtab[ihi]);
                             fb     = ftab[ihi]+slope*(b-xtab[ihi]);
                             result = 0.5f*(b-a)*(fa+fb);
                          }
                          else {
                             result = 0.0f;

                             for(i = ilo; i != ihi-1; ++i) {
                                 result = result+0.5f*(xtab[i+1]-xtab[i])
                                          * (ftab[i]+ftab[i+1]);
                             }

                             if(ilo==1) {
                                 slope  = (ftab[1]-ftab[0])/(xtab[1]-xtab[0]);
                                 fa     = ftab[0]+slope*(a-xtab[0]);
                                 result = result+0.5f*(xtab[ilo]-a)*(fa+ftab[ilo]);
                             }
                             else {
                                 slope  = (ftab[ilo]-ftab[ilo-1])/(xtab[ilo]-xtab[ilo-1]);
                                 fa     = ftab[ilo-1]+slope*(a-xtab[ilo-1]);
                                 result = result+0.5f*(xtab[ilo]-a)*(fa+ftab[ilo]);
                             }
 
                            if(ihi==ntab) {
                                 slope  = (ftab[ntab]-ftab[ntab-1])/(xtab[ntab]-xtab[ntab-1] ;
                                 fb     = ftab[ntab-1]+slope*(b-xtab[ntab-1]);
                                 result = result+0.5*(b-xtab[ntab])*(fb+ftab[ntab]);
                            }
                            else {
                                 slope  = (ftab[ihi+1]-ftab[ihi])/(xtab[ihi+1]-xtab[ihi]);
                                 fb     = ftab[ihi]+slope*(b-xtab[ihi]);
                                 result = result+0.5f*(b-xtab[ihi])*(fb+ftab[ihi]);
                            }
                         }

                         if(ind!=1) {
                             ind = 1;
                             syl = b;
                             b = a;
                             a = syl;
                             result = -result;
                          }
                          
                     } 
















 
      } // math

} // gms














#endif /*__GMS_PLINT_QUAD_HPP__*/
