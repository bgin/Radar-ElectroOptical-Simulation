

#ifndef __GMS_RCS_COMMON_ZMM16R4_HPP__
#define __GMS_RCS_COMMON_ZMM16R4_HPP__ 040120231012

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

    const unsigned int GMS_RCS_COMMON_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMMON_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMMON_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMMON_ZMM16R4_FULLVER =
      1000U*GMS_RCS_COMMON_ZMM16R4_MAJOR+
      100U*GMS_RCS_COMMON_ZMM16R4_MINOR+
      10U*GMS_RCS_COMMON_ZMM16R4_MICRO;
    const char * const GMS_RCS_COMMON_ZMM16R4_CREATION_DATE = "04-01-2023 10:12 AM +00200 (WED 04 01 2023 GMT+2)";
    const char * const GMS_RCS_COMMON_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMMON_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMMON_ZMM16R4_DESCRIPTION   = "AVX512 optimized Radar Cross Section common functions."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_zmm16r4.hpp"


namespace  gms {

        namespace radiolocation {


                   /*
                        Complex wavenumber, vector of 16 varying values of
                        complex permeability and permitivity.
                    */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void k_zmm16r4(const __m512 mur,
                                  const __m512 mui
                                  const __m512 epsr,
                                  const __m512 epsi,
                                  const __m512 om,
                                  __m512 * __restrict kr,
                                  __m512 * __restrict ki) {

                        using namespace gms::math;
                        __m512 sqrr,sqri;
                        __m512 t0r,t0i;
                        __m512 wrkc = _mm512_setzero_ps();
                        cmul_zmm16r4(mur,mui,epsr,epsi,&t0r,&t0i);
                        csqrt_zmm16r4(t0r,t0i,&wrkc,&sqrr,&sqri);
                        *kr = _mm512_mul_ps(om,sqrr);
                        *ki = _mm512_mul_ps(om,sqri);
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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rc_zmm16r4(const __m512 x,
                                     const __m512 y,
                                     const __m512 errtot) {

                          const register __m512 c1 = _mm512_set1_ps(0.142857142857142857142857142857f);
                          const register __m512 c2 = _mm512_set1_ps(0.409090909090909090909090909091f);
                          const register __m512 _1 = _mm512_set1_ps(1.0f);
                          const register __m512 thr= _mm512_set1_ps(0.333333333333333333333333333333f);
                          const register __m512 c3 = _mm512_set1_ps(0.375f);
                          const register __m512 _2 = _mm512_set1_ps(2.0f);
                          const register __m512 qtr= _mm512_set1_ps(0.25f);
                          const register __m512 c4 = _mm512_set1_ps(0.3f);
                          register __m512 rc,xn,yn,mu,s,sn,lamda;
                          xn = x;
                          yn = y;
                          
                          while(true) {

                                mu = _mm512_mul_ps(_mm512_add_ps(xn,
                                                          _mm512_add_ps(yn,yn)),thr);
                                sn = _mm512_div_ps(_mm512_add_ps(yn,mu),
                                                   _mm512_sub_ps(mu,_2));
                                if(_mm512_cmp_mask_ps(
                                             _mm512_abs_ps(sn),errtot,_CMP_LT_OQ)) {
                                   register __m512 t0 = _mm512_fmadd_ps(sn,c2,c3);
                                   register __m512 t1 = _mm512_fmadd_ps(t0,sn,c1);
                                   register __m512 t2 = _mm512_fmadd_ps(t1,sn,c4);
                                   s                  = _mm512_mul_ps(_mm512_mul_ps(sn,sn),t2);
                                   rc                 = _mm512_div_ps(_mm512_add_ps(_1,s),
                                                                      _mm512_sqrt_ps(mu));
                                   return (rc);
                                }
                                register __m512 sxn = _mm512_sqrt_ps(xn);
                                register __m512 syn = _mm512_sqrt_ps(yn);
                                lamda =  _mm512_fmadd_ps(_mm512_mul_ps(_2,sxn),syn,yn);
                                xn    = _mm512_mul_ps(_mm512_add_ps(xn,lamda),qtr);
                                yn    = _mm512_mul_ps(_mm512_add_ps(ym,lamda),qtr);
                         }
                }


                   

     } // radiolocation

} // gms













#endif /*__GMS_RCS_COMMON_ZMM16R4_HPP__*/

