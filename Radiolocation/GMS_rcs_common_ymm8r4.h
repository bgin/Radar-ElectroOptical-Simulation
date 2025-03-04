

#ifndef __GMS_RCS_COMMON_YMM8R4_H__
#define __GMS_RCS_COMMON_YMM8R4_H__ 210820240805

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_version {

    const unsigned int GMS_RCS_COMMON_YMM8R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMMON_YMM8R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMMON_YMM8R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMMON_YMM8R4_FULLVER =
      1000U*GMS_RCS_COMMON_YMM8R4_MAJOR+
      100U*GMS_RCS_COMMON_YMM8R4_MINOR+
      10U*GMS_RCS_COMMON_YMM8R4_MICRO;
    const char * const GMS_RCS_COMMON_YMM8R4_CREATION_DATE = "21-08-2024 08:05 PM +00200 (TUE 21 AUG 2024 GMT+2)";
    const char * const GMS_RCS_COMMON_YMM8R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMMON_YMM8R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMMON_YMM8R4_DESCRIPTION   = "AVX (single precision) optimized Radar Cross Section common functions."

}

#if !defined(__AVX512F__) || !defined(__AVX52VL__)
#error "Support of AVX512F or AVX512VL ISA required!!"
#endif

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_ymm8r4.hpp"



namespace  gms {

        namespace radiolocation {


                   /*
                        Complex wavenumber, vector of 16 varying values of
                        complex permeability and permitivity.
                    */
                   
	          
	          
                  
	           static inline
                   void k_ymm8r4(const __m256 mur,
                                  const __m256 mui
                                  const __m256 epsr,
                                  const __m256 epsi,
                                  const __m256 om,
                                  __m256 * __restrict kr,
                                  __m256 * __restrict ki) {

                        using namespace gms::math;
                        __m256 sqrr,sqri;
                        __m256 t0r,t0i;
                        __m256 wrkc = _mm256_setzero_ps();
                        cmul_ymm8c4(mur,mui,epsr,epsi,&t0r,&t0i);
                        csqrt_ymm8c4(t0r,t0i,&wrkc,&sqrr,&sqri);
                        *kr = _mm256_mul_ps(om,sqrr);
                        *ki = _mm256_mul_ps(om,sqri);
               }


                 /*

c*********************************************************************72
c
cc RC computes the elementary integral RC(X,Y).
c
c  Discussion:
c
c    This function computes the elementary integral
c
c      RC(X,Y) = Integral ( 0 <= T < oo )
c
c                              -1/2     -1
c                    (1/2)(T+X)    (T+Y)  DT,
c
c    where X is nonnegative and Y is positive.  The duplication
c    theorem is iterated until the variables are nearly equal,
c    and the function is then expanded in Taylor series to fifth
c    order.  
c
c    Logarithmic, inverse circular, and inverse hyperbolic 
c    functions can be expressed in terms of RC.  
c
c    Check by addition theorem: 
c
c      RC(X,X+Z) + RC(Y,Y+Z) = RC(0,Z),
c      where X, Y, and Z are positive and X * Y = Z * Z.
c
c  Modified:
c
c    27 May 2018
c
c  Author:
c
c    Bille Carlson, Elaine Notis
c
c  Reference:
c
c    Bille Carlson,
c    Computing Elliptic Integrals by Duplication,
c    Numerische Mathematik,
c    Volume 33, 1979, pages 1-16.
c
c    Bille Carlson, Elaine Notis,
c    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 3, pages 398-403, September 1981.
c
c  Parameters:
c
c    Input, double precision X, Y, the arguments in the integral.
c
c    Input, double precision ERRTOL, the error tolerance.
c    Relative error due to truncation is less than
c      16 * ERRTOL ^ 6 / (1 - 2 * ERRTOL).
c    Sample choices:  
c      ERRTOL   Relative truncation error less than
c      1.D-3    2.D-17
c      3.D-3    2.D-14
c      1.D-2    2.D-11
c      3.D-2    2.D-8
c      1.D-1    2.D-5
c
c    Output, integer IERR, the error flag.
c    0, no error occurred.
c    1, abnormal termination.
c
                    */


                   
	          
	          
                  
	            __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 rc_ymm8r4(const __m256 x,
                                     const __m256 y,
                                     const __m256 errtot); 


                   
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 rc_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) px,
                                       const float * __restrict __ATTR_ALIGN__(32) py,
                                       const float * __restrict __ATTR_ALIGN__(32) perrtot); 


                    
	          
	          
                  
	              __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 rc_ymm8r4_u(const float * __restrict  px,
                                       const float * __restrict  py,
                                       const float * __restrict  perrtot); 


                       


                 /*

c*********************************************************************72
c
cc RD computes an incomplete elliptic integral of the second kind, RD(X,Y,Z).
c
c  Discussion:
c
c    This function computes an incomplete elliptic integral of the second kind.
c
c    RD(X,Y,Z) = Integral ( 0 <= T < oo )
c
c                                -1/2     -1/2     -3/2
c                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,
c
c    where X and Y are nonnegative, X + Y is positive, and Z is positive.
c
c    If X or Y is zero, the integral is complete.
c
c    The duplication theorem is iterated until the variables are
c    nearly equal, and the function is then expanded in Taylor
c    series to fifth order.  
c
c    Check: 
c
c      RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) = 3 / sqrt ( X * Y * Z ), 
c      where X, Y, and Z are positive.
c
c  Modified:
c
c    27 May 2018
c
c  Author:
c
c    Bille Carlson, Elaine Notis
c
c  Reference:
c
c    Bille Carlson,
c    Computing Elliptic Integrals by Duplication,
c    Numerische Mathematik,
c    Volume 33, 1979, pages 1-16.
c
c    Bille Carlson, Elaine Notis,
c    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 3, pages 398-403, September 1981.
c
c  Parameters:
c
c    Input, double precision X, Y, Z, the arguments in the integral.
c
c    Input, double precision ERRTOL, the error tolerance.
c    The relative error due to truncation is less than
c      3 * ERRTOL ^ 6 / (1-ERRTOL) ^ 3/2.
c    Sample choices:
c      ERRTOL   Relative truncation error less than
c      1.D-3    4.D-18
c      3.D-3    3.D-15
c      1.D-2    4.D-12
c      3.D-2    3.D-9
c      1.D-1    4.D-6
c
c    Output, integer IERR, the error flag.
c    0, no error occurred.
c    1, abnormal termination.
c
                   */  


                   
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 rd_ymm8r4(const __m256 x,
                                     const __m256 y,
                                     const __m256 z,
                                     const __m256 errtot);

                    
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 rd_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) px,
                                       const float * __restrict __ATTR_ALIGN__(32) py,
                                       const float * __restrict __ATTR_ALIGN__(32) pz,
                                       const float * __restrict __ATTR_ALIGN__(32) perrtot);

                      
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 rd_ymm8r4_u(const float * __restrict  px,
                                       const float * __restrict  py,
                                       const float * __restrict  pz,
                                       const float * __restrict  perrtot); 

/*							fresnl.c
 *
 *	Fresnel integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, S, C;
 * void fresnl();
 *
 * fresnl( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
 *
 * The integrals are evaluated by a power series for x < 1.
 * For x >= 1 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
 * S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
 *
 *
 *
 * ACCURACY:
 *
 *  Relative error.
 *
 * Arithmetic  function   domain     # trials      peak         rms
 *   IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
 *   IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
 *   DEC        S(x)      0, 10        6000       2.2e-16     3.9e-17
 *   DEC        C(x)      0, 10        5000       2.3e-16     3.9e-17
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/





                   
	          
	          
                  
	                 __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   void fresnel_ymm8r4(const __m256 xxa,
                                        __m256 * __restrict ssa,
                                        __m256 * __restrict cca); 

                   
	          
	          
                  
	                __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   void fresnel_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pxxa,
                                          float * __restrict __ATTR_ALIGN__(32) ssa,
                                          float * __restrict __ATTR_ALIGN__(32) cca); 

                   
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   void fresnel_ymm8r4_u(const float * __restrict  pxxa,
                                          float * __restrict  ssa,
                                          float * __restrict  cca); 



                  /*
                           Same as above -- divided into Fresnel 'C' integral
                           and Fresnel 'S' integral.
                     */


                   
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 fresnel_C_ymm8r4(const __m256 xxa); 

                   
	          
	          
                  
	              __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 fresnel_C_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pxxa);


                  
	          
	          
                  
	              __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 fresnel_C_ymm8r4_u(const float * __restrict  pxxa); 

                   
	          
	          
                       __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 fresnel_S_ymm8r4(const __m256 xxa); 


                    
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 fresnel_S_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pxxa); 


                   
	          
	          
                  
	               __ATTR_HOT__
	            __ATTR_VECTORCALL__
                   __m256 fresnel_S_ymm8r4_u(const float * __restrict  pxxa); 

                   





     } // radiolocation

} // gms













#endif /*__GMS_RCS_COMMON_YMM8R4_H__*/

