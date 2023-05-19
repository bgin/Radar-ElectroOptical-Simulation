

#ifndef __GMS_SPEC_FUN_ZMM8R8_HPP__
#define __GMS_SPEC_FUN_ZMM8R8_HPP__


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

    const unsigned int GMS_SPEC_FUN_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_SPEC_FUN_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_SPEC_FUN_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_SPEC_FUN_ZMM8R8_FULLVER =
      1000U*GMS_SPEC_FUN_ZMM8R8_MAJOR+
      100U*GMS_SPEC_FUN_ZMM8R8_MINOR+
      10U*GMS_SPEC_FUN_ZMM8R8_MICRO;
    const char * const GMS_SPEC_FUNC_ZMM8R8_CREATION_DATE = "18-05-2023 08:38 AM +00200 (THR 18 05 2023 GMT+2)";
    const char * const GMS_SPEC_FUNC_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_SPEC_FUNC_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_SPEC_FUNC_ZMM8R8_DESCRIPTION   = "Special Functions library manually vectorized (avx512 double)."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.hpp"
#include "GMS_sleefsimddp.hpp"




namespace gms {


        namespace math {
        
        
/*
               !*****************************************************************************80
!
!! BESEI0 evaluates the exponentially scaled Bessel I0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the modified Bessel
!    function of the first kind of order zero multiplied by EXP(-ABS(X)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESEI0, the value of the function.
!
               
*/
        
        
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besei0_zmm8r8(const __m512d x) {
	           
	                  __m512d result;
	                  int32_t jint;
	                  jint = 2;
	                  result = calci0_zmm8r8(x,jint);
	                  return (result); 
	         }  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besei0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                  register __m512d x = _mm512_load_pd(&px[0]);
	                  register __m512d result;
	                  int32_t jint;
	                  jint = 2;
	                  result = calci0_zmm8r8(x,jint);
	                  return (result); 
	         }  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besei0_zmm8r8_u(const double * __restrict  px) {
	           
	                  register __m512d x = _mm512_loadu_pd(&px[0]);
	                  register __m512d result;
	                  int32_t jint;
	                  jint = 2;
	                  result = calci0_zmm8r8(x,jint);
	                  return (result); 
	         }  
	         
	         
/*
!*****************************************************************************80
!
!! BESEI1 evaluates the exponentially scaled Bessel I1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the first kind of order one
!    multiplied by EXP(-ABS(X)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESEI1, the value of the function.
*/	         


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besei1_zmm8r8(const __m512 x) {
	           
	                  register __m512d result;
	                  int32_t jint;
	                  jint = 2;
	                  result = calci1_zmm8r8(x,jint);
	                  return (result); 
	         }  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besei1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                  register __m512d x = _mm512_load_pd(&px[0]);
	                  register __m512d result;
	                  int32_t jint;
	                  jint = 2;
	                  result = calci1_zmm8r8(x,jint);
	                  return (result); 
	         }  
	         
	         
	         
        
       } // math




} // gms































#endif /*__GMS_SPEC_FUNC_ZMM8R8_HPP__*/
