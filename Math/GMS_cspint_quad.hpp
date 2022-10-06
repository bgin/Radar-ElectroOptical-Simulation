

#ifndef __GMS_CSPINT_QUAD_HPP__
#define __GMS_CSPINT_QUAD_HPP__ 051020221039


namespace file_info {


        const unsigned int GMS_CSPINT_QUAD_MAJOR = 1;

	const unsigned int GMS_CSPINT_QUAD_MINOR = 1;

	const unsigned int GMS_CSPINT_QUAD_MICRO = 0;

	const unsigned int GMS_CSPINT_QUAD_FULLVER = 
		1000U*GMS_CSPINT_QUAD_MAJOR+100U*GMS_CSPINT_QUAD_MINOR+10U*GMS_CSPINT_QUAD_MICRO;

	const char * const GMS_CSPINT_QUAD_CREATE_DATE = "05-10-2022 10:36 +00200 (WED 05 OCT 2022 GMT+2)";

	const char * const GMS_CSPINT_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_CSPINT_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}

/*
!*****************************************************************************80
!
!! CSPINT estimates the integral of a tabulated function.
!
!  Discussion:
!
!    The routine is given the value of a function F(X) at a set of 
!    nodes XTAB, and estimates
!
!      Integral ( A <= X <= B ) F(X) DX
!
!    by computing the cubic natural spline S(X) that interpolates
!    F(X) at the nodes, and then computing
!
!      Integral ( A <= X <= B ) S(X) DX
!
!    exactly.
!
!    Other output from the program includes the definite integral
!    from X(1) to X(I) of S(X), and the coefficients necessary for
!    the user to evaluate the spline S(X) at any point.
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
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), contains the points at which the
!    function was evaluated.  The XTAB's must be distinct and
!    in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated values of
!    the function, FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, lower limit of integration.
!
!    Input, real ( kind = 8 ) B, upper limit of integration.
!
!    Output, real ( kind = 8 ) Y(3,NTAB), will contain the coefficients
!    of the interpolating natural spline over each subinterval.
!    For XTAB(I) <= X <= XTAB(I+1),
!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!                   + Y(2,I)*(X-XTAB(I))**2
!                   + Y(3,I)*(X-XTAB(I))**3
!
!    Output, real ( kind = 8 ) E(NTAB), E(I) = the definite integral from
!    XTAB(1) to XTAB(I) of S(X).
!
!    Workspace, real ( kind = 8 ) WORK(NTAB).
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
*/

#include <omp.h>
#include <cstdint>
#include <limits> //Nan value
#include "GMS_config.h"

namespace gms {

      namespace math {
          
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void cspint(  const int32_t ntab,
                                      const double * __restrict __ATTR_ALIGN__(64) xtab
                                      const double * __restrict __ATTR_ALIGN__(64) ftab
                                      const double A,
                                      const double B,
                                      double * __restrict __ATTR_ALIGN__(64) y1,
                                      double * __restrict __ATTR_ALIGN__(64) y2,
                                      double * __restrict __ATTR_ALIGN__(64) y3,
                                      double * __restrict __ATTR_ALIGN__(64) e,
                                      double * __restrict __ATTR_ALIGN__(64) work,
                                      double & result) {
                           
                           if(ntab<3) {
                               result = std::numeric_limits<double>::quiet_NaN();
                               return;
                           }
                           double r,s,term,u,t0;
                           int32_t i,j;

                           __assume_aligned(xtab,64);
                           __assume_aligned(ftab,64);
                           __assume_aligned(y1,64);
                           __assume_aligned(y2,64);
                           __assume_aligned(y3,64);
                           __assume_aligned(e,64);
                           __assume_aligned(work,64);

                           for(i = 0; i != ntab-1; ++i){
                               if(xtab[i+1]<=xtab[i]) {
                                  result = std::numeric_limits<double>::quiet_NaN();
                                  return;
                               }
                           }
                           s = 0.0;
                           #pragma vector aligned
                           #pragma omp simd
                           for(i = 0; i != ntab-1; ++i) {
                               r     = (ftab[i+1]-ftab[i])/(xtab[i+1]-xtab[i]);
                               y2[i] = r-s;
                               s     = r;
                           }
                           s = 0.0;
                           r = 0.0;
                           y2[0]    = 0.0;
                           y2[ntab] = 0.0;

                           for(i = 1; i != ntab-1; ++i) {
                               y2[i]  = y2[i]+r*y2[i-1];
                               work[i]= 2.0*(xtab[i-1]-xtab[i+1])-r*s;
                               s      = xtab[i+1]-xtab[i];
                               r      = s/work[i];
                           }

                           #pragma vector aligned
                           #pragma omp simd
                           for(j = 1; j != ntab-1; ++j) {
                               i = ntab+1-j;
                               y2[i] = (xtab[i+1]-xtab[i]*y2[i+1]-y2[i])/work[i];
                           }

                           #pragma vector aligned
                           #pragma omp simd
                           for(i = 0; i != ntab-1; ++i) {
                               s     = xtab[i+1]-xtab[i];
                               r     = y2[i+1]-y2[i];
                               y3[i] = r/s;
                               y2[i] = 3.0*y2[i];
                               t0    = y2[i]+r;
                               y1[i] = (ftab[i+1]-ftab[i])/(s-t0*s);
                           }

                           e[0] = 0.0;
                           #pragma vector aligned
                           #pragma omp simd reduction(+:e)
                           for(i = 0; i != ntab-1; ++i) {
                               s    = xtab[i+1]-xtab[i];
                               term = ((( y3[i]*0.25*s+y2[i]*0.33333333333333333333333333333)*s
                                      + y1[i]*0.5)*s+ftab[i])*s;
                               e[i+1] = e[i]+term; 
                           }
                          //!
                          //!  Determine where the endpoints A and B lie in the mesh of XTAB's.
                          //!
                          r = A;
                          u = 1.0;
                          result = 0.0;
                          for(j = 0, j != 1; ++j) {
                              //The endpoint is less than or equal to XTAB(1).
                              if(r<=xtab[0]) {
                                 result = result-u*((r-xtab[0])* y1[0]*0.5+
                                 ftab[0])*(r-xtab[0]); 
                              }
                              else if(xtab[ntab]<=r) {
                              // The endpoint is greater than or equal to XTAB(NTAB).
                                result = result-u*(e[ntab]+(r-xtab[ntab]) 
                                         *(ftab[ntab]+0.5*(ftab[ntab-1] 
                                         +(xtab[ntab]-xtab[ntab-1])*y1[ntab-1]) 
                                         *(r-xtab[ntab])));
                              }
                              else {
                               //  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
                                 for(i = 0; i != ntab-1; ++i) {
                                     if(r<=xtab[i+1]) {
                                        r = r-xtab[i];
                                        result = result-u*(e[i]+(((
                                                 y3[i]*0.25*r 
                                                 +y2[i]*0.333333333333333333333333333333333)*r
                                                 +y1[i]*0.5)*r+ftab[i])*r);
                                        goto label120;
                                     }
                                 }

                              }
label120:
                                 u = -1.0;
                                 r = B;
                          }
                      }


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void cspint(  const int32_t ntab,
                                      const float * __restrict __ATTR_ALIGN__(64) xtab
                                      const float * __restrict __ATTR_ALIGN__(64) ftab
                                      const float A,
                                      const float B,
                                      float * __restrict __ATTR_ALIGN__(64) y1,
                                      float * __restrict __ATTR_ALIGN__(64) y2,
                                      float * __restrict __ATTR_ALIGN__(64) y3,
                                      float * __restrict __ATTR_ALIGN__(64) e,
                                      float * __restrict __ATTR_ALIGN__(64) work,
                                      float & result) {
                           
                           if(ntab<3) {
                              result = std::numeric_limits<float>::quiet_NaN();
                              return;
                           }
                           float r,s,term,u,t0;
                           int32_t i,j;

                           __assume_aligned(xtab,64);
                           __assume_aligned(ftab,64);
                           __assume_aligned(y1,64);
                           __assume_aligned(y2,64);
                           __assume_aligned(y3,64);
                           __assume_aligned(e,64);
                           __assume_aligned(work,64);

                           for(i = 0; i != ntab-1; ++i) {
                               if(xtab[i+1]<=xtab[i]) {
                                  result = std::numeric_limits<float>::quiet_NaN();
                                  return;
                               }
                           }
                           s = 0.0f;
                           #pragma vector aligned
                           #pragma omp simd
                           for(i = 0; i != ntab-1; ++i) {
                               r     = (ftab[i+1]-ftab[i])/(xtab[i+1]-xtab[i]);
                               y2[i] = r-s;
                               s     = r;
                           }
                           s = 0.0f;
                           r = 0.0f;
                           y2[0]    = 0.0f;
                           y2[ntab] = 0.0f;

                           for(i = 1; i != ntab-1; ++i) {
                               y2[i]  = y2[i]+r*y2[i-1];
                               work[i]= 2.0f*(xtab[i-1]-xtab[i+1])-r*s;
                               s      = xtab[i+1]-xtab[i];
                               r      = s/work[i];
                           }

                           #pragma vector aligned
                           #pragma omp simd
                           for(j = 1; j != ntab-1; ++j) {
                               i = ntab+1-j;
                               y2[i] = (xtab[i+1]-xtab[i]*y2[i+1]-y2[i])/work[i];
                           }

                           #pragma vector aligned
                           #pragma omp simd
                           for(i = 0; i != ntab-1; ++i) {
                               s     = xtab[i+1]-xtab[i];
                               r     = y2[i+1]-y2[i];
                               y3[i] = r/s;
                               y2[i] = 3.0f*y2[i];
                               t0    = y2[i]+r;
                               y1[i] = (ftab[i+1]-ftab[i])/(s-t0*s);
                           }

                           e[0] = 0.0f;
                           #pragma vector aligned
                           #pragma omp simd reduction(+:e)
                           for(i = 0; i != ntab-1; ++i) {
                               s    = xtab[i+1]-xtab[i];
                               term = ((( y3[i]*0.25f*s+y2[i]*0.33333333333333333333333333333f)*s
                                      + y1[i]*0.5f)*s+ftab[i])*s;
                               e[i+1] = e[i]+term; 
                           }
                          //!
                          //!  Determine where the endpoints A and B lie in the mesh of XTAB's.
                          //!
                          r = A;
                          u = 1.0f;
                          result = 0.0f;
                          for(j = 0, j != 1; ++j) {
                              //The endpoint is less than or equal to XTAB(1).
                              if(r<=xtab[0]) {
                                 result = result-u*((r-xtab[0])* y1[0]*0.5f+
                                 ftab[0])*(r-xtab[0]); 
                              }
                              else if(xtab[ntab]<=r) {
                              // The endpoint is greater than or equal to XTAB(NTAB).
                                result = result-u*(e[ntab]+(r-xtab[ntab]) 
                                         *(ftab[ntab]+0.5f*(ftab[ntab-1] 
                                         +(xtab[ntab]-xtab[ntab-1])*y1[ntab-1]) 
                                         *(r-xtab[ntab])));
                              }
                              else {
                               //  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
                                 for(i = 0; i != ntab-1; ++i) {
                                     if(r<=xtab[i+1]) {
                                        r = r-xtab[i];
                                        result = result-u*(e[i]+(((
                                                 y3[i]*0.25f*r 
                                                 +y2[i]*0.333333333333333333333333333333333f)*r
                                                 +y1[i]*0.5f)*r+ftab[i])*r);
                                        goto label120;
                                     }
                                 }

                              }
label120:
                                 u = -1.0f;
                                 r = B;
                          }
                      }




                      


      } //math


}// gms







#endif /*__GMS_CSPINT_QUAD_HPP__*/
