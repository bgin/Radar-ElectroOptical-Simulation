

#ifndef __GMS_FILON_SIN_QUAD_HPP__
#define __GMS_FILON_SIN_QUAD_HPP__ 071020220903


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
                              __ATTR_HOT__
                              __ATTR_ALIGN__(32)
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
                              __ATTR_HOT__
                              __ATTR_ALIGN__(32)
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

#include "GMS_malloc.h"


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void filon_sin(const int32_t ntab,
                                       const double * __restrict __ATTR_ALIGN__(64) ftab,
                                       const double a,
                                       const double b,
                                       const double t,
                                       double &result) {
                           using namespace gms::common;
                           if(__builtin_expect(a==b,0)) {
                              result = std::numeric_limits<double>::quiet_NaN();
                              return;
                           }
                           if(__builtin_expect(ntab<=1,0) || 
                              __builtin_expect((ntab%2)!=1,0)) {
                               result = std::numeric_limits<double>::quiet_NaN();
                               return;
                            }
                            double alpha,beta,c2n,c2nm1;
                            double s2n,s2nm1,ftabn,xtabn;
                            double cost,gamma,h,sint;
                            double * __restrict xtab = NULL;
                            double theta,th2,th3,th4,th5,th6;
                            double th7,th8,ftab1,xtab1;
                            double sum,t0,t1,t2,t3,t4;
                            const double ntab1 = (double)(ntab-1);
                            int32_t i;
                            xtab = (double*)gms_mm_malloc((std::size_t)ntab,64ULL);
                            vec_even(ntab,a,b,xtab);
                            ftabn= ftab[ntab];
                            ftab1= ftab[0];
                            xtabn= xtab[ntab];
                            xtab1= xtab[0];
                            h    = (b-a)/ntab1;
                            theta= t*h;
                            th2  = theta*theta;
                            th3  = th2*theta;
                            sint = std::sin(theta);
                            th4  = th3*theta;
                            th5  = th4*theta;
                            cost = std::cos(theta);
                            th6  = th5*theta;
                            th7  = th6*theta;
                            th8  = th7*theta;
                            t0   = std::sin(t*xtabn);
                            t1   = std::sin(t*xtab1);
                            if(6.0*std::abs(theta)<=1.0) {
#if (DIVISION_REPLACEMENT) == 0
                                 alpha = 2.0*th3/45.0-2.0*th5/315.0+
                                         2.0*th7/4725.0;
                                 beta  = 0.666666666666666666666666666667+
                                         2.0*th2/15.0-4.0*th4/105.0+
                                         2.0*th6/567.0-4.0*th8/22275.0;
                                 gamma = 1.333333333333333333333333333333-
                                         2.0*th2/15.0+th4/210.0-th6/11340.0;
         
                    
#else
                                 alpha = 2.0*th3*0.022222222222222222222222222222222-
                                         2.0*th5*0.002857142857142857142857142857+
                                         2.0*th7*0.000211640211640211640211640212;
                                 beta  = 2.0*0.333333333333333333333333333333333+
                                         2.0*th2*0.066666666666666666666666666667-
                                         4.0*th4*0.009523809523809523809523809524+
                                         2.0*th6*0.00176366843033509700176366843-
                                         4.0*th8*0.000044893378226711560044893378;
                                 gamma = 4.0*0.333333333333333333333333333333333333-
                                         2.0*th2*0.066666666666666666666666666667+
                                         th4*0.004761904761904761904761904762-
                                         th6*0.000088183421516754850088183422; 
#endif
                            }
                            else {
                                
                                alpha = (th2+theta*sint*cost-2.0*sint*sint)/th3;
                                beta  = (2.0*theta+2.0*theta*cost*cost-4.0*sint*cost)/th3;
                                gamma = 4.0*(sint-theta*cost)/th3;
                             }
                            t2    = std::sin(t*xtab1);
                            t3    = std::sin(t*xtabn);
                            s2n   = 0.0;
                            s2nm1 = 0.0;
                            sum   = 0.0;
                            __assume_aligned(ftab,64);
                            __assume_aligned(xtab,64);
                            #pragma vector aligned
                            #pragma omp simd reduction(+:sum)
                            for(i = 0; i != ntab; i += 2) {
                                const double ftabi = ftab[i];
                                const double xtabi = xtab[i];
                                sum = sum+ftabi*std::sin(t*xtabi);
                                s2n = sum-0.5*(ftabn*t0+ftab1*t1);
                            }
                            #pragma vector aligned
                            #pragma omp simd reduction(+:s2nm1)
                            for(i = 1; i != ntab-1; i += 2) {
                                const double ftabi = ftab[i];
                                const double xtabi = xtab[i];
                                s2nm1 = s2nm1+ftabi*std::sin(t*xtabi);
                            }
                            sum    = beta*s2n+gamma*s2nm1;
                            t4     = ftab1*t2-ftab2*t3;
                            result = alpha*t4+sum;
                            gms_mm_free(xtab);
                      }  


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			static inline
                        void filon_sin(const int32_t ntab,
                                       const float * __restrict __ATTR_ALIGN__(64) ftab,
                                       const float a,
                                       const float b,
                                       const float t,
                                       double &result) {
                           using namespace gms::common;
                           if(__builtin_expect(a==b,0)) {
                              result = std::numeric_limits<float>::quiet_NaN();
                              return;
                           }
                           if(__builtin_expect(ntab<=1,0) || 
                              __builtin_expect((ntab%2)!=1,0)) {
                               result = std::numeric_limits<float>::quiet_NaN();
                               return;
                            }
                            float alpha,beta;
                            float s2n,s2nm1,ftabn,xtabn;
                            float cost,gamma,h,sint;
                            float * __restrict xtab = NULL;
                            float theta,th2,th3,th4,th5,th6;
                            float th7,th8,ftab1,xtab1;
                            float sum,t0,t1,t2,t3,t4;
                            const float ntab1 = (float)(ntab-1);
                            int32_t i;
                            xtab = (float*)gms_mm_malloc((std::size_t)ntab,64ULL);
                            vec_even(ntab,a,b,xtab);
                            ftabn= ftab[ntab];
                            ftab1= ftab[0];
                            xtabn= xtab[ntab];
                            xtab1= xtab[0];
                            h    = (b-a)/ntab1;
                            theta= t*h;
                            th2  = theta*theta;
                            th3  = th2*theta;
                            sint = std::sin(theta);
                            th4  = th3*theta;
                            th5  = th4*theta;
                            cost = std::cos(theta);
                            th6  = th5*theta;
                            th7  = th6*theta;
                            th8  = th7*theta;
                            t0   = std::sin(t*xtabn);
                            t1   = std::sin(t*xtab1);
                            if(6.0f*std::abs(theta)<=1.0f) {
#if (DIVISION_REPLACEMENT) == 0
                                 alpha = 2.0f*th3/45.0f-2.0f*th5/315.0f+
                                         2.0f*th7/4725.0;
                                 beta  = 0.666666666666666666666666666667f+
                                         2.0f*th2/15.0f-4.0f*th4/105.0f+
                                         2.0f*th6/567.0f-4.0f*th8/22275.0f;
                                 gamma = 1.333333333333333333333333333333f-
                                         2.0f*th2/15.0f+th4/210.0f-th6/11340.0f;
         
                    
#else
                                 alpha = 2.0f*th3*0.022222222222222222222222222222222f-
                                         2.0f*th5*0.002857142857142857142857142857f+
                                         2.0f*th7*0.000211640211640211640211640212f;
                                 beta  = 2.0f*0.333333333333333333333333333333333f+
                                         2.0f*th2*0.066666666666666666666666666667f-
                                         4.0f*th4*0.009523809523809523809523809524f+
                                         2.0f*th6*0.00176366843033509700176366843f-
                                         4.0f*th8*0.000044893378226711560044893378f;
                                 gamma = 4.0f*0.333333333333333333333333333333333333f-
                                         2.f0*th2*0.066666666666666666666666666667f+
                                         th4*0.004761904761904761904761904762f-
                                         th6*0.000088183421516754850088183422f; 
#endif
                            }
                            else {
                                
                                alpha = (th2+theta*sint*cost-2.0f*sint*sint)/th3;
                                beta  = (2.0f*theta+2.0f*theta*cost*cost-4.0f*sint*cost)/th3;
                                gamma = 4.0*(sint-theta*cost)/th3;
                             }
                            t2    = std::sin(t*xtab1);
                            t3    = std::sin(t*xtabn);
                            s2n   = 0.0f;
                            s2nm1 = 0.0f;
                            sum   = 0.0f;
                            __assume_aligned(ftab,64);
                            __assume_aligned(xtab,64);
                            #pragma vector aligned
                            #pragma omp simd reduction(+:sum)
                            for(i = 0; i != ntab; i += 2) {
                                const double ftabi = ftab[i];
                                const double xtabi = xtab[i];
                                sum = sum+ftabi*std::sin(t*xtabi);
                                s2n = sum-0.5f*(ftabn*t0+ftab1*t1);
                            }
                            #pragma vector aligned
                            #pragma omp simd reduction(+:s2nm1)
                            for(i = 1; i != ntab-1; i += 2) {
                                const double ftabi = ftab[i];
                                const double xtabi = xtab[i];
                                s2nm1 = s2nm1+ftabi*std::sin(t*xtabi);
                            }
                            sum    = beta*s2n+gamma*s2nm1;
                            t4     = ftab1*t2-ftab2*t3;
                            result = alpha*t4+sum;
                            gms_mm_free(xtab);
                      }  

                                  
            

     } // math

} // gms






#endif /*__GMS_FILON_SIN_QUAD_HPP__*/
