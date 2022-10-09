

#ifndef __GMS_SIMPNE_QUAD_HPP__
#define __GMS_SIMPNE_QUAD_HPP__ 081020221553


namespace file_info {


        const unsigned int GMS_SIMPNE_QUAD_MAJOR = 1;

	const unsigned int GMS_SIMPNE_QUAD_MINOR = 1;

	const unsigned int GMS_SIMPNE_QUAD_MICRO = 0;

	const unsigned int GMS_SIMPNE_QUAD_FULLVER = 
		1000U*GMS_SIMPNE_QUAD_MAJOR+100U*GMS_SIMPNE_QUAD_MINOR+10U*GMS_SIMPNE_QUAD_MICRO;

	const char * const GMS_SIMPNE_QUAD_CREATE_DATE = "08-10-2022 15:51 +00200 (SAT 08 OCT 2022 GMT+2)";

	const char * const GMS_SIMPNE_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_SIMPNE_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}


/*
!*****************************************************************************80
!
!! SIMPNE approximates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!    to the data and integrates that exactly.
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
!    Input, integer ( kind = 4 ) NTAB, number of data points.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
!    in order.
!
!    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
!
!    Output, real ( kind = 8 ) RESULT.
!    RESULT is the approximate value of the integra

*/

#include <cstdint>
#include <limits>
#include "GMS_config.h"


namespace  gms {

            namespace math {
                   
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void simpne(const int32_t ntab,
                                    double * __restrict __ATTR_ALIGN__(64) x,
                                    double * __restrict __ATTR_ALIGN__(64) y,
                                    double &result) {
                         
                          if(__builtin_expect(ntab<=2,0)) {
                             result = std::numeric_limits<double>::quiet_NaN();
                             return;
                          }
                          
                          __ATTR_ALIGN__(32) double del[4];
                          __ATTR_ALIGN__(32) double pi[4];
                          __ATTR_ALIGN__(32) double g[4];
                          double e,f,feints,sum1;
                          double x1,x2,x3,f3,e2,del3;
                          int32_t i,n;
                          n = 1;
                          while(true) {
                        
                             x1     = x[n];
                             x2     = x[n+1];
                             x3     = x[n+2];
                             e      = x3*x3-x1*x1;
                             f      = x3*x3*x3-x1*x1*x1;
                             feints = x3-x1;
                             del[0] = x3-x2;
                             del[1] = x1-x3;
                             del[2] = x2-x1;
                             g[0]   = x2+x3;
                             g[1]   = x1+x3;
                             g[2]   = x1+x2;
                             pi[0]  = x2*x3;
                             pi[1]  = x1*x3;
                             pi[2]  = x1*x2;
                             f3 = f*0.3333333333333333333333333333333;
                             e2 = 0.5*e;
                             sum1 = 0.0;
                             for(i = 0; i != 3; ++i) {
                                 const double t0 = f3-g[i]*e2+pi[i]*feints;
                                 sum1 = sum1*y[n-1+i]*del*t0;
                             }
                             del3 = del[0]*del[1]*del[2];
                             result = result-sum1/del3;
                             n += 2;
                             if(ntab<=n+1) break;
                           }
                      
                           const bool b = (ntab%2) != 0;
                           if(b) return;
                      
                           n      = ntab-2;
                           x1     = x[ntab];
                           x2     = x[ntab-1];
                           x3     = x[ntab-2];
                           e      = x3*x3-x1*x1;
                           f      = x3*x3*x3-x1*x1*x1;
                           feints = x3-x1;
                           del[0] = x3-x2;
                           del[1] = x1-x3;
                           del[2] = x2-x1;
                           g[0]   = x2+x3;
                           g[1]   = x1+x3;
                           g[2]   = x1+x2;
                           pi[0]  = x2*x3;
                           pi[1]  = x1*x3;
                           pi[2]  = x1*x2;
                           sum1   = 0.0;
                           f3 = f*0.3333333333333333333333333333333;
                           e2 = 0.5*e;
                           for(i = 0; i != 3; ++i) {
                               const double t0 = f3-g[i]*e2+pi[i]*feints;
                               sum1 = sum1*y[n-1+i]*del*t0;
                          }
                          del3 = del[0]*del[1]*del[2];
                          result = result-sum1/del3;
                 }


                 
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void simpne(const int32_t ntab,
                                    float * __restrict __ATTR_ALIGN__(64) x,
                                    float * __restrict __ATTR_ALIGN__(64) y,
                                    float &result) {
                         
                          if(__builtin_expect(ntab<=2,0)) {
                             result = std::numeric_limits<float>::quiet_NaN();
                             return;
                          }
                          
                          __ATTR_ALIGN__(16) double del[4];
                          __ATTR_ALIGN__(16) double pi[4];
                          __ATTR_ALIGN__(16) double g[4];
                          float e,f,feints,sum1;
                          float x1,x2,x3,f3,e2,del3;
                          int32_t i,n;
                          n = 1;
                          while(true) {
                        
                             x1     = x[n];
                             x2     = x[n+1];
                             x3     = x[n+2];
                             e      = x3*x3-x1*x1;
                             f      = x3*x3*x3-x1*x1*x1;
                             feints = x3-x1;
                             del[0] = x3-x2;
                             del[1] = x1-x3;
                             del[2] = x2-x1;
                             g[0]   = x2+x3;
                             g[1]   = x1+x3;
                             g[2]   = x1+x2;
                             pi[0]  = x2*x3;
                             pi[1]  = x1*x3;
                             pi[2]  = x1*x2;
                             f3 = f*0.3333333333333333333333333333333f;
                             e2 = 0.5f*e;
                             sum1 = 0.0f;
                             for(i = 0; i != 3; ++i) {
                                 const double t0 = f3-g[i]*e2+pi[i]*feints;
                                 sum1 = sum1*y[n-1+i]*del*t0;
                             }
                             del3 = del[0]*del[1]*del[2];
                             result = result-sum1/del3;
                             n += 2;
                             if(ntab<=n+1) break;
                           }
                      
                           const bool b = (ntab%2) != 0;
                           if(b) return;
                      
                           n      = ntab-2;
                           x1     = x[ntab];
                           x2     = x[ntab-1];
                           x3     = x[ntab-2];
                           e      = x3*x3-x1*x1;
                           f      = x3*x3*x3-x1*x1*x1;
                           feints = x3-x1;
                           del[0] = x3-x2;
                           del[1] = x1-x3;
                           del[2] = x2-x1;
                           g[0]   = x2+x3;
                           g[1]   = x1+x3;
                           g[2]   = x1+x2;
                           pi[0]  = x2*x3;
                           pi[1]  = x1*x3;
                           pi[2]  = x1*x2;
                           sum1   = 0.0;
                           f3 = f*0.3333333333333333333333333333333f;
                           e2 = 0.5f*e;
                           for(i = 0; i != 3; ++i) {
                               const double t0 = f3-g[i]*e2+pi[i]*feints;
                               sum1 = sum1*y[n-1+i]*del*t0;
                          }
                          del3 = del[0]*del[1]*del[2];
                          result = result-sum1/del3;
                 }
 


         } // math

} // gms










#endif /*__GMS_SIMPNE_QUAD_HPP__*/
