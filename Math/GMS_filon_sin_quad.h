

#ifndef __GMS_FILON_SIN_QUAD_H__
#define __GMS_FILON_SIN_QUAD_H__ 071020220903


namespace file_info {


        const unsigned int GMS_FILON_SIN_QUAD_MAJOR = 1;

	const unsigned int GMS_FILON_SIN_QUAD_MINOR = 1;

	const unsigned int GMS_FILON_SIN_QUAD_MICRO = 0;

	const unsigned int GMS_FILON_SIN_QUAD_FULLVER = 
		1000U*GMS_FILON_SIN_QUAD_MAJOR+100U*GMS_FILON_SIN_QUAD_MINOR+10U*GMS_FILON_SIN_QUAD_MICRO;

	const char * const GMS_FILON_SIN_QUAD_CREATE_DATE = "08-10-2022 10:16 +00200 (SAT 08 OCT 2022 GMT+2)";

	const char * const GMS_FILON_SIN_QUAD_BUILD_DATE = __DATE__ " " __TIME__;

	const char * const GMS_FILON_SIN_QUAD_AUTHOR = "Tabulated data integrator ported from Fortran version";

}


/*
  *****************************************************************************80
!
!! FILON_COS uses Filon's method on integrals with a cosine factor.
!
!  Discussion:
!
!    The integral to be approximated has the form:
!
!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!
!    where T is user specified.
!
!    The function is interpolated over each subinterval by
!    a parabolic arc.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Chase, Lloyd Fosdick,
!    An Algorithm for Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, Number 8, August 1969, pages 453-457.
!
!    Stephen Chase, Lloyd Fosdick,
!    Algorithm 353:
!    Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, Number 8, August 1969, pages 457-458.
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
!    Input, integer ( kind = 4 ) NTAB, the number of data points.
!    NTAB must be odd, and greater than 1.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the value of the function
!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) T, the multiplier of the X argument of the cosine.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!    

*/


#include <cstdint>
#include <limits> //Nan value
#include <cmath>
#include "GMS_config.h"


namespace  gms {

            namespace  math {


                         namespace {

                              __ATTR_ALWAYS_INLINE__
                            
			      static inline
                              void vec_even(const int32_t n,
                                            const double ahi,
                                            const double alo,
                                            double * __restrict __ATTR_ALIGN__(64) a) {
                                    
                                    if(__builtin_expect(n==1,0)) {
                                       a[0] = 0.5*(alo+ahi);
                                    }
                                    else {
                                       __assume_aligned(a,64);
                                       const double n1 = (double)(n-1);
                                       for(int32_t i = 1; i != n; ++i) {
                                           const double ni = (double)(n-i);
                                           const double i1 = (double)(i-1);
                                           a[i] = ni*alo+i1*ahi/n1;
                                       }
                                    }
                               }

                              __ATTR_ALWAYS_INLINE__
                             
			      static inline
                              void vec_even(const int32_t n,
                                            const float ahi,
                                            const float alo,
                                            float * __restrict __ATTR_ALIGN__(64) a) {
                                    
                                    if(__builtin_expect(n==1,0)) {
                                       a[0] = 0.5f*(alo+ahi);
                                    }
                                    else {
                                       __assume_aligned(a,64);
                                       const float n1 = (double)(n-1);
                                       for(int32_t i = 1; i != n; ++i) {
                                           const double ni = (double)(n-i);
                                           const double i1 = (double)(i-1);
                                           a[i] = ni*alo+i1*ahi/n1;
                                       }
                                    }
                               }
                               
                         }




                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
                        void filon_sin(const int32_t ntab,
                                       const double * __restrict __ATTR_ALIGN__(64) ftab,
                                       const double a,
                                       const double b,
                                       const double t,
                                       double &result); 
                                       

                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		
                        void filon_sin(const int32_t ntab,
                                       const float * __restrict __ATTR_ALIGN__(64) ftab,
                                       const float a,
                                       const float b,
                                       const float t,
                                       double &result); 
                                  
            

     } // math

} // gms






#endif /*__GMS_FILON_SIN_QUAD_H__*/
