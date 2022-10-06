
#ifndef __GMS_CUBINT_QUAD_HPP__
#define __GMS_CUBINT_QUAD_HPP__ 061020220904


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
#include <limits> //Nan value
#include <algorithm>
#include "GMS_config.h"

namespace  gms {

        namespace  math {

                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void cubint(const int32_t ntab,
                                    double * __restrict __ATTR_ALIGN__(64) xtab,
                                    double * __restrict __ATTR_ALIGN__(64) ftab,
                                    const int32_t ia,
                                    const int32_t ib,
                                    double &result,
                                    double &error) {
                         
                          if(__builtin_expect(ia==ib,0)) {
                            result = std::numeric_limits<double>::quiet_NaN();
                            error  = result;
                            return;
                         }
                         if(__builtin_expect(ntab<4,0) || 
                            __builtin_expect(ntab<ib,0)) {
                            result = std::numeric_limits<double>::quiet_NaN();
                            error  = result;
                            return;
                         }
                         if(__builtin_expect(ia<1,0) || 
                            __builtin_expect(ib<1,0)) {
                            result = std::numeric_limits<double>::quiet_NaN();
                            error  = result;
                            return;
                         }
                         if(__builtin_expect(ntab<ia,0)) {
                            result = std::numeric_limits<double>::quiet_NaN();
                            error  = result;
                            return;
                         }
                        
                         double c,d1,d2,d3;
                         double h1,h2,h3,h4;
                         double r1,r2,r3,r4;
                         double s,term;
                         int32_t i,ia,ib,ind;
                         int32_t it,j,k;
                         bool b0;
                         /*
                              !  Temporarily switch IA and IB, and store minus sign in IND
                              !  so that, while integration is carried out from low X's
                              !  to high ones, the sense of the integral is preserved.
                          */
                          if(ib<ia) {
                             ind = -1;
                             it = ib;
                             ib = ia
                             ia = it
                          }
                          else {
                             ind = 1;
                          }
                          s = 0.0;
                          c = 0.0;
                          r4= 0.0;
                          j = ntab-2;
                          b0 = ntab==4;
                          if(ia<ntab-1 || b0) {
                             j = std::max(3,ia);
                          }
                          k = 4;
                          if(2<ib || b0) {
                             k = std::min(ntab,ib+2)-1;
                          } 
 
                          for(i = j; i != k; ++i) {
                              if(i<=j) {
                                 h2 = xtab[j-1]-xtab[j-2];
                                 d3 = (ftab[j-1]-ftab[j-2])/h2;
                                 h3 = xtab[j]-xtab[j-1];
                                 d1 = (ftab[j]-ftab[j-1])/h3;
                                 h1 = h2 + h3;
                                 d2 = (d1-d3)/h1;
                                 h4 = xtab[j+1]-xtab[j];
                                 r1 = (ftab[j+1]-ftab[j])/h4;
                                 r2 = (r1-d1)/(h4+h3);
                                 h1 = h1 + h4;
                                 r3 = (r2-d2)/h1;
                                 if(ia<=1) {
                                    result = h2*(ftab[1]+h2*(0.5*d3-h2 
                                             *(d2/6.0-(h2+h3+h3)*r3/12.0)));
                                              s = -h2*h2*h2*(h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0;
                                 }  
                              }
                              else {
                                  h4 = xtab[i+1]-xtab[i];
                                  r1 = (ftab[i+1]-ftab[i])/h4;
                                  r4 = h4+h3;
                                  r2 = (r1-d1)/r4;
                                  r4 = r4+h2;
                                  r3 = (r2-d2)/r4;
                                  r4 = (r3-d3)/(r4+h1);
                              }
                              
                              if(ia<i && i<=ib) {
                                 term = h3*((ftab[i]+ftab[i-1])*0.5 
                                       -h3*h3*(d2+r2+(h2-h4)*r3)/12.0 ;
                                 result = result+term;
                                 c = h3*h3*h3*(2.0*h3*h3
                                     + 5.0*(h3*(h4+h2)+2.0*h2*h4))/120.0;
                                 error = error+(c+s)*r4;
                                 
                                 if(i!=j)
                                    s = c;
                                 else
                                    s = s+c+c;
                              }
                              else {
                                
                                 error = error+r4*s;
                              }
                             
                              if(k<=i) {

                                 if(ntab<=ib) {
                                    term = h4*(ftab[ntab]-h4*(0.5*r1
                                           + h4*(r2/6.0+(h3+h3+h4)*r3/12.0)));
                                    result = result+term;
                                    error  = error-h4*h4*h4*r4*
                                             (h4*(3.0*h4+5.0*h2) 
                                             + 10.0*h3*(h2+h3+h4))/60.0;
                                    
                                 }
                                 
                                 if(ntab-1<=ib) error = error+s*r4;
                              }
                              else {
                                
                                 h1 = h2;
                                 h2 = h3;
                                 h3 = h4;
                                 d1 = r1;
                                 d2 = r2;
                                 d3 = r3;
                              }
                              
                          } 
                          /*
                               !  Restore original values of IA and IB, reverse signs
                               !  of RESULT and ERROR, to account for integration
                               !  that proceeded from high X to low X. 
                           */
                          if(ind!=1) {
                              it = ib;
                              ib = ia;
                              ia = it;
                              result = -result;
                              error = -error;
                          }
                  }


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void cubint(const int32_t ntab,
                                    float * __restrict __ATTR_ALIGN__(64) xtab,
                                    float * __restrict __ATTR_ALIGN__(64) ftab,
                                    const int32_t ia,
                                    const int32_t ib,
                                    float &result,
                                    float &error) {
                         
                          if(__builtin_expect(ia==ib,0)) {
                            result = std::numeric_limits<float>::quiet_NaN();
                            error  = result;
                            return;
                         }
                         if(__builtin_expect(ntab<4,0) || 
                            __builtin_expect(ntab<ib,0)) {
                            result = std::numeric_limits<float>::quiet_NaN();
                            error  = result;
                            return;
                         }
                         if(__builtin_expect(ia<1,0) || 
                            __builtin_expect(ib<1,0)) {
                            result = std::numeric_limits<float>::quiet_NaN();
                            error  = result;
                            return;
                         }
                         if(__builtin_expect(ntab<ia,0)) {
                            result = std::numeric_limits<float>::quiet_NaN();
                            error  = result;
                            return;
                         }
                        
                         float c,d1,d2,d3;
                         float h1,h2,h3,h4;
                         float r1,r2,r3,r4;
                         float s,term;
                         int32_t i,ia,ib,ind;
                         int32_t it,j,k;
                         bool b0;
                         /*
                              !  Temporarily switch IA and IB, and store minus sign in IND
                              !  so that, while integration is carried out from low X's
                              !  to high ones, the sense of the integral is preserved.
                          */
                          if(ib<ia) {
                             ind = -1;
                             it = ib;
                             ib = ia
                             ia = it
                          }
                          else {
                             ind = 1;
                          }
                          s = 0.0f;
                          c = 0.0f;
                          r4= 0.0f;
                          j = ntab-2;
                          b0 = ntab==4;
                          if(ia<ntab-1 || b0) {
                             j = std::max(3,ia);
                          }
                          k = 4;
                          if(2<ib || b0) {
                             k = std::min(ntab,ib+2)-1;
                          } 
 
                          for(i = j; i != k; ++i) {
                              if(i<=j) {
                                 h2 = xtab[j-1]-xtab[j-2];
                                 d3 = (ftab[j-1]-ftab[j-2])/h2;
                                 h3 = xtab[j]-xtab[j-1];
                                 d1 = (ftab[j]-ftab[j-1])/h3;
                                 h1 = h2 + h3;
                                 d2 = (d1-d3)/h1;
                                 h4 = xtab[j+1]-xtab[j];
                                 r1 = (ftab[j+1]-ftab[j])/h4;
                                 r2 = (r1-d1)/(h4+h3);
                                 h1 = h1 + h4;
                                 r3 = (r2-d2)/h1;
                                 if(ia<=1) {
                                    result = h2*(ftab[1]+h2*(0.5f*d3-h2 
                                             *(d2/6.0f-(h2+h3+h3)*r3/12.0f)));
                                              s = -h2*h2*h2*(h2*(3.0f*h2+5.0f*h4)+10.0f*h3*h1)/60.0f;
                                 }  
                              }
                              else {
                                  h4 = xtab[i+1]-xtab[i];
                                  r1 = (ftab[i+1]-ftab[i])/h4;
                                  r4 = h4+h3;
                                  r2 = (r1-d1)/r4;
                                  r4 = r4+h2;
                                  r3 = (r2-d2)/r4;
                                  r4 = (r3-d3)/(r4+h1);
                              }
                              
                              if(ia<i && i<=ib) {
                                 term = h3*((ftab[i]+ftab[i-1])*0.5f 
                                       -h3*h3*(d2+r2+(h2-h4)*r3)/12.0f ;
                                 result = result+term;
                                 c = h3*h3*h3*(2.0f*h3*h3
                                     + 5.0f*(h3*(h4+h2)+2.0*h2*h4))/120.0f;
                                 error = error+(c+s)*r4;
                                 
                                 if(i!=j)
                                    s = c;
                                 else
                                    s = s+c+c;
                              }
                              else {
                                
                                 error = error+r4*s;
                              }
                             
                              if(k<=i) {

                                 if(ntab<=ib) {
                                    term = h4*(ftab[ntab]-h4*(0.5f*r1
                                           + h4*(r2/6.0f+(h3+h3+h4)*r3/12.0f)));
                                    result = result+term;
                                    error  = error-h4*h4*h4*r4*
                                             (h4*(3.0f*h4+5.0f*h2) 
                                             + 10.0f*h3*(h2+h3+h4))/60.0f;
                                    
                                 }
                                 
                                 if(ntab-1<=ib) error = error+s*r4;
                              }
                              else {
                                
                                 h1 = h2;
                                 h2 = h3;
                                 h3 = h4;
                                 d1 = r1;
                                 d2 = r2;
                                 d3 = r3;
                              }
                              
                          } 
                          /*
                               !  Restore original values of IA and IB, reverse signs
                               !  of RESULT and ERROR, to account for integration
                               !  that proceeded from high X to low X. 
                           */
                          if(ind!=1) {
                              it = ib;
                              ib = ia;
                              ia = it;
                              result = -result;
                              error = -error;
                          }
                  }

  
                  
       }
}






















#endif /*__GMS_CUBINT_QUAD_HPP__*/
