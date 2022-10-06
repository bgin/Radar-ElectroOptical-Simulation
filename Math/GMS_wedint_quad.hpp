

#ifndef __GMS_WEDINT_QUAD_HPP__
#define __GMS_WEDINT_QUAD_HPP__ 061020221649


namespace file_info {


        const unsigned int GMS_WEDINT_QUAD_MAJOR = 1;

	const unsigned int GMS_WEDINT_QUAD_MINOR = 1;

	const unsigned int GMS_WEDINT_QUAD_MICRO = 0;

	const unsigned int GMS_WEDINT_QUAD_FULLVER = 
		1000U*GMS_WEDINT_QUAD_MAJOR+100U*GMS_WEDINT_QUAD_MINOR+10U*GMS_WEDINT_QUAD_MICRO;

	const char * const GMS_WEDINT_QUAD_CREATE_DATE = "06-10-2022 16:49 +00200 (THR 06 OCT 2022 GMT+2)";

	const char * const GMS_WEDINT_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_WEDINT_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}

/*
 *****************************************************************************80
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, is the number of data points.  
!    (NTAB-1) must be divisible by 6.
!
!    Input, real ( kind = 8 ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated data values.
!
!    Output, real ( kind = 8 ) RESULT, is the approximation to the integral.   
*/

#include <cstdint>
#include <omp.h>
#include <limits> //Nan value
#include "GMS_config.h"


namespace  gms {

         namespace  math {
              
                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void wedint(const int32_t ntab,
                                    const double h,
                                    const double * __restrict __ATTR_ALIGN__(64) ftab,
                                    double &result) {

                           if(__builtin_expect(ntab<=1,0)) {
                                result = std::numeric_limits<double>::quiet_NaN();
                                return;
                             }
                           if((ntab%6)!=1) {
                                result = std::numeric_limits<double>::quiet_NaN();
                                return;
                           }
                          int32_t i;
                          __assume_aligned(ftab,64);
                          #pragma vector aligned
                          #pragma omp simd reduction(+:result)
                          for(i = 0; i != ntab-6; i += 6) {
                              result = result+ftab[i]+5.0*ftab[i+1]+
                                       ftab[i+2]+6.0*ftab[i+3]     +
                                       ftab[i+4]+5.0*ftab[i+5]     +
                                       ftab[i+6];
                          }
                          result = 3.0*h*result/10.0;
                     }


                       __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void wedint(const int32_t ntab,
                                    const float h,
                                    const float * __restrict __ATTR_ALIGN__(64) ftab,
                                    double &result) {

                           if(__builtin_expect(ntab<=1,0)) {
                                result = std::numeric_limits<float>::quiet_NaN();
                                return;
                             }
                           if((ntab%6)!=1) {
                                result = std::numeric_limits<float>::quiet_NaN();
                                return;
                           }
                          int32_t i;
                          __assume_aligned(ftab,64);
                          #pragma vector aligned
                          #pragma omp simd reduction(+:result)
                          for(i = 0; i != ntab-6; i += 6) {
                              result = result+ftab[i]+5.0f*ftab[i+1]+
                                       ftab[i+2]+6.0f*ftab[i+3]     +
                                       ftab[i+4]+5.0f*ftab[i+5]     +
                                       ftab[i+6];
                          }
                          result = 3.0f*h*result/10.0f;
                     }
     }

}










#endif /*__GMS_WEDINT_QUAD_HPP__*/
