

#ifndef __GMS_HIORDQ_QUAD_HPP__
#define __GMS_HIORDQ_QUAD_HPP__ 081020221251


namespace file_info {


        const unsigned int GMS_HIORDQ_QUAD_MAJOR = 1;

	const unsigned int GMS_HIORDQ_QUAD_MINOR = 1;

	const unsigned int GMS_HIORDQ_QUAD_MICRO = 0;

	const unsigned int GMS_HIORDQ_QUAD_FULLVER = 
		1000U*GMS_HIORDQ_QUAD_MAJOR+100U*GMS_HIORDQ_QUAD_MINOR+10U*GMS_HIORDQ_QUAD_MICRO;

	const char * const GMS_HIORDQ_QUAD_CREATE_DATE = "08-10-2022 12:51 +00200 (SAT 08 OCT 2022 GMT+2)";

	const char * const GMS_HIORDQ_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_HIORDQ_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}

/*
!*****************************************************************************80
!
!! HIORDQ approximates the integral of a function using equally spaced data.
!
!  Discussion:
!
!    The method applies the trapezoidal rule to various subsets of the
!    data, and then applies Richardson extrapolation.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    Alan Kaylor Cline,
!    Department of Computer Science,
!    University of Texas at Austin.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, number of data points.
!
!    Input, real ( kind = 8 ) DELT, the spacing between the X values of the
!    data.  The actual X values are not needed!
!
!    Input, real ( kind = 8 ) Y(NTAB), the Y values of the data.
!
!    Work array, real ( kind = 8 ) WORK(2*(NTAB-1)).  The actual minimum amount
!    of workspace required is two times the number of integer
!    divisors of NTAB-1.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.

*/

#include <cstdint>
#include "GMS_config.h"


namespace gms {

          namespace math {
                
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void hiordq(const int32_t ntab,
                                    const double delt,
                                    const double * __restrict __ATTR_ALIGN__(64) y,
                                    const double * __restrict __ATTR_ALIGN__(64) work,
                                    double &result) {

                            double fac,sum2,sum1;
                            int32_t i,j,jbak,jj,k;
                            int32_t ntab1,ntab1k;
                            ntab1 = ntab-1;
                            
                            sum1 = 0.5*(y[0]+y[ntab]);
                            j = -1;

                            for(k = 1; k != ntab1; ++k) {
                                ntab1k = ntab1/k;
                                if(ntab1k*k == ntab1) {
                                   sum = -sum1;
                                   for(i = 0; i != ntab; i += ntab1k) {
                                       sum2 += y[i];
                                   }
                                   j += 2;
                                   const double t0 = (double)ntab1k;
                                   work[j]   = delt*sum2*t0;
                                   work[j+1] = t0*t0;

                                   if(k != 1) {
                                      for(jj = 3; jj != j; jj += 2) {
                                          jbak = j+1-jj;
                                          const double wj1 = work[j+1];
                                          fac  = wj1/(wj1-work[jbak+1]);
                                          work[jbak] = work[jbak+2]+fac*(work[jbak]-work[jbak+2]);
                                      }
                                   }
                                }
                            }
                            result = work[0];
                    }


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void hiordq(const int32_t ntab,
                                    const float delt,
                                    const float * __restrict __ATTR_ALIGN__(64) y,
                                    const float * __restrict __ATTR_ALIGN__(64) work,
                                    float &result) {

                            float fac,sum2,sum1;
                            int32_t i,j,jbak,jj,k;
                            int32_t ntab1,ntab1k;
                            ntab1 = ntab-1;
                            
                            sum1 = 0.5f*(y[0]+y[ntab]);
                            j = -1;

                            for(k = 1; k != ntab1; ++k) {
                                ntab1k = ntab1/k;
                                if(ntab1k*k == ntab1) {
                                   sum = -sum1;
                                   for(i = 0; i != ntab; i += ntab1k) {
                                       sum2 += y[i];
                                   }
                                   j += 2;
                                   const float t0 = (float)ntab1k;
                                   work[j]   = delt*sum2*t0;
                                   work[j+1] = t0*t0;

                                   if(k != 1) {
                                      for(jj = 3; jj != j; jj += 2) {
                                          jbak = j+1-jj;
                                          const float wj1 = work[j+1];
                                          fac  = wj1/(wj1-work[jbak+1]);
                                          work[jbak] = work[jbak+2]+fac*(work[jbak]-work[jbak+2]);
                                      }
                                   }
                                }
                            }
                            result = work[0];
                    }


   
        } // math

} //gms





#endif /*__GMS_HIORDQ_QUAD_HPP__*/
