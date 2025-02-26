


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


#include <omp.h>
#include <limits> //Nan value
#include "GMS_wedint_quad.h"



                        void gms::math::wedint(const int32_t ntab,
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


                      
                        void gms::math::wedint(const int32_t ntab,
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
  
