
#ifndef __GMS_SIMPSN_QUAD_HPP__
#define __GMS_SIMPSN_QUAD_HPP__ 091020221053


namespace file_info {


        const unsigned int GMS_SIMPN_QUAD_MAJOR = 1;

	const unsigned int GMS_SIMPN_QUAD_MINOR = 1;

	const unsigned int GMS_SIMPN_QUAD_MICRO = 0;

	const unsigned int GMS_SIMPN_QUAD_FULLVER = 
		1000U*GMS_SIMPN_QUAD_MAJOR+100U*GMS_SIMPN_QUAD_MINOR+10U*GMS_SIMPN_QUAD_MICRO;

	const char * const GMS_SIMPN_QUAD_CREATE_DATE = "09-10-2022 10:53 +00200 (SUN 09 OCT 2022 GMT+2)";

	const char * const GMS_SIMPN_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_SIMPN_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}


/*

!*****************************************************************************80
!
!! SIMPSN approximates the integral of evenly spaced data.
!
!  Discussion:
!
!    Simpson's rule is used.
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
!    Input, integer ( kind = 4 ) NTAB, the number of data points.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) H, specifies the increment between the
!    X values.  Note that the actual X values are not needed,
!    just the constant spacing!
!
!    Input, real ( kind = 8 ) Y(NTAB), the data.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral
!    from the first to the last point.
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
                        void simpn(const int32_t ntab,
                                   const double h,
                                   double * __restrict __ATTR_ALIGN__(64) y,
                                   double &result) {

                            if(__builtin_expect(ntab<=2,0)) {
                             result = std::numeric_limits<double>::quiet_NaN();
                             return;
                           }
                          __ATTR_ALIGN__(32) double del[4];
                          __ATTR_ALIGN__(32) double pi[4];
                          __ATTR_ALIGN__(32) double g[4];
                          double f,e,sum1,f3,h2;
                          int32_t i,n;
                          bool b;
                          
                          b = (ntab%2)==0;
                          if(b)
                            n = ntab-1;
                          else
                            n = ntab;

                          result = y[0]+y[n]+4.0*y[n-1];
                          for(i = 1; i != n-1; i += 2) {
                              const double yi = y[i];
                              const double yi1= y[i+1];
                              result = result+4.0*yi+2.0*yi1;
                          }
                          result = h*result/3.0;
                          f      = h*h*h;
                          del[0] = h;
                          del[1] = -2.0*h;
                          del[2] = h;
                          g[0]   = h;
                          g[1]   = 0.0;
                          g[2]   = -h;
                          pi[0]  = 0.0;
                          pi[1]  = -h*h;
                          pi[2]  = 0.0;
                          n      = n-1;
                          sum1 = 0.0;
                          f3   = f*0.3333333333333333333333333333333;
                          h2   = 0.5*h*h;
                          for(i = 0; i != 2; ++i) {
                              const double t0 = f3-g[i]*h2+pi[i]*h;
                              sum1 = sum1+y[n-1+i]*del[i]*t0;
                          }
                          result = result+0.5*sum1/f;
                     }


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void simpn(const int32_t ntab,
                                   const double h,
                                   float * __restrict __ATTR_ALIGN__(64) y,
                                   float &result) {

                            if(__builtin_expect(ntab<=2,0)) {
                             result = std::numeric_limits<float>::quiet_NaN();
                             return;
                           }
                          __ATTR_ALIGN__(16) double del[4];
                          __ATTR_ALIGN__(16) double pi[4];
                          __ATTR_ALIGN__(16) double g[4];
                          float f,e,sum1,f3,h2;
                          int32_t i,n;
                          bool b;
                          
                          b = (ntab%2)==0;
                          if(b)
                            n = ntab-1;
                          else
                            n = ntab;

                          result = y[0]+y[n]+4.0f*y[n-1];
                          for(i = 1; i != n-1; i += 2) {
                              const double yi = y[i];
                              const double yi1= y[i+1];
                              result = result+4.0f*yi+2.0f*yi1;
                          }
                          result = h*result/3.0f;
                          f      = h*h*h;
                          del[0] = h;
                          del[1] = -2.0f*h;
                          del[2] = h;
                          g[0]   = h;
                          g[1]   = 0.0f;
                          g[2]   = -h;
                          pi[0]  = 0.0f;
                          pi[1]  = -h*h;
                          pi[2]  = 0.0f;
                          n      = n-1;
                          sum1 = 0.0f;
                          f3   = f*0.3333333333333333333333333333333f;
                          h2   = 0.5f*h*h;
                          for(i = 0; i != 2; ++i) {
                              const double t0 = f3-g[i]*h2+pi[i]*h;
                              sum1 = sum1+y[n-1+i]*del[i]*t0;
                          }
                          result = result+0.5f*sum1/f;
                     }


    } // math

} // gms


#endif /*__GMS_SIMPSN_QUAD_HPP__*/
