

#ifndef __GMS_CSPINT_QUAD_H__
#define __GMS_CSPINT_QUAD_H__ 051020221039


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


#include <cstdint>
#include "GMS_config.h"

namespace gms {

      namespace math {
      
                         /*
                   Work (input) arrays for integrator 'cspint'
                   Multi-threaded i.e. (two-threaded)
               */
             
                  __ATTR_ALIGN__(64) struct CSPINT_DATA_R8_2T {
               
                       double * __restrict  Ya1; 
                       double * __restrict  Ya2; 
                       double * __restrict  Ya3; 
                       double * __restrict  Ea;  
                       double * __restrict  WRKa; 
                       double * __restrict  Yb1;
                       double * __restrict  Yb2; 
                       double * __restrict  Yb3; 
                       double * __restrict  Eb; 
                       double * __restrict  WRKb;  
               };
               
                /*
                   Work (input) arrays for integrator 'cspint'
                   Single-threaded.
               */
                __ATTR_ALIGN__(64) struct CSPINT_DATA_R8_1T {
               
                       double * __restrict  Y1; 
                       double * __restrict  Y2; 
                       double * __restrict  Y3; 
                       double * __restrict  E;  
                       double * __restrict  WRK; 
                       
               };
          
                      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
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
                                      double & result); 
                                      

               __ATTR_ALIGN__(64) struct CSPINT_DATA_R4_2T {
               
                       float * __restrict  Ya1; 
                       float * __restrict  Ya2; 
                       float * __restrict  Ya3; 
                       float * __restrict  Ea;  
                       float * __restrict  WRKa; 
                       float * __restrict  Yb1;
                       float * __restrict  Yb2; 
                       float * __restrict  Yb3; 
                       float * __restrict  Eb; 
                       float * __restrict  WRKb;  
               };
               
                /*
                   Work (input) arrays for integrator 'cspint'
                   Single-threaded.
               */
                __ATTR_ALIGN__(64) struct CSPINT_DATA_R4_1T {
               
                       float * __restrict  Y1; 
                       float * __restrict  Y2; 
                       float * __restrict  Y3; 
                       float * __restrict  E;  
                       float * __restrict  WRK; 
                       
               }; 

                       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
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
                                      float & result); 



                      


      } //math


}// gms







#endif /*__GMS_CSPINT_QUAD_H__*/
