

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
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besei1_zmm8r8_u(const double * __restrict  px) {
	           
	                  register __m512d x = _mm512_loadu_pd(&px[0]);
	                  register __m512d result;
	                  int32_t jint;
	                  jint = 2;
	                  result = calci1_zmm8r8(x,jint);
	                  return (result); 
	         }  
	         
	     
/*
   !*****************************************************************************80
!
!! BESEK0 evaluates the exponentially scaled Bessel K0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    multiplied by the exponential function.
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
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!    0 < X.
!
!    Output, real ( kind = 8 ) BESK0, the value of the function.
*/


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besek0_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 2;
	                   result = calck0_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besek0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 2;
	                   result = calck0_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besek0_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 2;
	                   result = calck0_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
/*
!*****************************************************************************80
!
!! BESEK1 evaluates the exponentially scaled Bessel K1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order one
!    multiplied by the exponential function, for arguments
!    XLEAST <= ARG <= XMAX.
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
!    Output, real ( kind = 8 ) BESEK1, the value of the function.
*/	  


                   
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besek1_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 2;
	                   result = calck1_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besek1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 2;
	                   result = calck1_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besek1_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 2;
	                   result = calck1_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	         
/*
     !*****************************************************************************80
!
!! BESI0 evaluates the Bessel I0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for
!    modified Bessel functions of the first kind of order zero for
!    arguments ABS(ARG) <= XMAX.
!
!    See comments heading CALCI0.
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
!    Output, real ( kind = 8 ) BESI0, the value of the function.
*/


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi0_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calci0_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calci0_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi0_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512 x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calci0_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
/*
    !*****************************************************************************80
!
!! BESI1 evaluates the Bessel I1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for
!    modified Bessel functions of the first kind of order one for
!    arguments ABS(ARG) <= XMAX.
!
!    See comments heading CALCI1.
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
!    Output, real ( kind = 8 ) BESI1, the value of the function.        
*/


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi1_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calci1_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calci1_zmm8r8(x,jint);
	                   return (result);
	          }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi1_zmm8r8_u(const double * __restrict px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calci1_zmm8r8(x,jint);
	                   return (result);
	          }
	          
	          
/*
    *****************************************************************************80
!
!! BESJ0 evaluates the Bessel J0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the first kind of order zero for arguments  |X| <= XMAX
!
!    See comments heading CALJY0.
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
!    Output, real ( kind = 8 ) BESJ0, the value of the function.       
*/


                
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besj0_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 0;
	                   result = caljy0_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besj0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 0;
	                   result = caljy0_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	              __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besj0_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 0;
	                   result = caljy0_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
/*
*****************************************************************************80
!
!! BESJ1 evaluates the Bessel J1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the first kind of order zero for arguments  |X| <= XMAX
!
!    See comments heading CALJY1.
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
!    Output, real ( kind = 8 ) BESJ1, the value of the function.
*/	 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besj1_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 0;
	                   result = caljy1_zmm8r8(x,jint);
	                   return (result);
	          }     
	          
	          
	             
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besj1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 0;
	                   result = caljy1_zmm8r8(x,jint);
	                   return (result);
	          }    
	          
	          
	             __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besj1_zmm8r8_u(const double * __restrict px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 0;
	                   result = caljy1_zmm8r8(x,jint);
	                   return (result);
	          }  
	          
	          

/*
    !*****************************************************************************80
!
!! BESK0 evaluates the Bessel K0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    for arguments 0.0 < ARG <= XMAX.
!
!    See comments heading CALCK0.
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
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) BESK0, the value of the function.
*/	


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besk0_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calck0_zmm8r8(x,jint);
	                   return (result);
	          }     
	          
	          
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besk0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calck0_zmm8r8(x,jint);
	                   return (result);
	          }     
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besk0_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calck0_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	          
/*
      *****************************************************************************80
!
!! BESK1 evaluates the Bessel K1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order one
!    for arguments XLEAST <= ARG <= XMAX.
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
!    Output, real ( kind = 8 ) BESK1, the value of the function.        
*/  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besk1_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calck1_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besk1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calck1_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besk1_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = calck1_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
/*
  !*****************************************************************************80
!
!! BESY0 evaluates the Bessel Y0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the second kind of order zero for arguments 0 < X <= XMAX.
!
!    See comments heading CALJY0.
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
!    Output, real ( kind = 8 ) BESY0, the value of the function.
*/	      


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besy0_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = caljy0_zmm8r8(x,jint);
	                   return (result);
	          }     
	          
	          
	            __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besy0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = caljy0_zmm8r8(x,jint);
	                   return (result);
	          }    
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besy0_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = caljy0_zmm8r8(x,jint);
	                   return (result);
	          }    
	          
	          
/*
    !*****************************************************************************80
!
!! BESY1 evaluates the Bessel Y1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for Bessel functions
!    of the second kind of order zero for arguments 0 < X <= XMAX.
!
!    See comments heading CALJY1.
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
!    Output, real ( kind = 8 ) BESY1, the value of the function.     

*/    


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besy1_zmm8r8(const __m512d x) {
	           
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = caljy1_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besy1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px) {
	           
	                   register __m512d x = _mm512_load_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = caljy1_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besy1_zmm8r8_u(const double * __restrict  px) {
	           
	                   register __m512d x = _mm512_loadu_pd(&px[0]);
	                   register __m512d result;
	                   int32_t jint;
	                   jint = 1;
	                   result = caljy1_zmm8r8(x,jint);
	                   return (result);
	          }   
	          
	         
/*
         !*****************************************************************************80
!
!! CALCEI computes various exponential integrals.
!
!  Discussion:
!
!    This routine computes the exponential integrals Ei(x),
!    E1(x), and  exp(-x)*Ei(x) for real arguments x where
!
!           integral (from t=-oo to t=x) (exp(t)/t),  x > 0,
!    Ei(x) =
!          -integral (from t=-x to t=+oo) (exp(t)/t),  x < 0,
!
!    and where the first integral is a principal value integral.
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
!  Reference:
!
!    William Cody, Henry Thacher,
!    Rational Chebyshev Approximations for the Exponential
!    Integral E1(x),
!    Mathematics of Computation,
!    Volume 22, Number 103, July 1968, pages 641-649.
!
!    William Cody, Henry Thacher,
!    Chebyshev Approximations for the Exponential
!    Integral Ei(x),
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 289-303.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  The argument must not
!    be zero.  If JINT = 2, then the argument must be strictly positive.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = EI ( ARG );
!    2, RESULT = EONE ( ARG );
!    3, RESULT = EXPEI ( ARG ).
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, Ei(x);
!    2, -Ei(-x);
!    3, exp(-x)*Ei(x).
!
*/


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calcei_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                                 
	                   __ATTR_ALIGN__(64) static const __m512d a[7]  = {_mm512_set1_pd(1.1669552669734461083368e+2),
	                                                       _mm512_set1_pd(2.1500672908092918123209e+3), 
                                                               _mm512_set1_pd(1.5924175980637303639884e+4), 
                                                               _mm512_set1_pd(8.9904972007457256553251e+4), 
                                                               _mm512_set1_pd(1.5026059476436982420737e+5),
                                                               _mm512_set1_pd(-1.4815102102575750838086e+5), 
                                                               _mm512_set1_pd(5.0196785185439843791020e+0)};
                           __ATTR_ALIGN__(64) static const __m512d b[6]  = {_mm512_set1_pd(4.0205465640027706061433e+1), 
                                                               _mm512_set1_pd(7.5043163907103936624165e+2),
                                                               _mm512_set1_pd(8.1258035174768735759855e+3), 
                                                               _mm512_set1_pd(5.2440529172056355429883e+4), 
                                                               _mm512_set1_pd(1.8434070063353677359298e+5), 
                                                               _mm512_set1_pd(2.5666493484897117319268e+5)};
                           __ATTR_ALIGN__(64) static const __m512d c[9]  = {_mm512_set1_pd(3.828573121022477169108e-1), 
                                                               _mm512_set1_pd(1.107326627786831743809e+1), 
                                                               _mm512_set1_pd(7.246689782858597021199e+1), 
                                                               _mm512_set1_pd(1.700632978311516129328e+2), 
                                                               _mm512_set1_pd(1.698106763764238382705e+2), 
                                                               _mm512_set1_pd(7.633628843705946890896e+1), 
                                                               _mm512_set1_pd(1.487967702840464066613e+1), 
                                                               _mm512_set1_pd(9.999989642347613068437e-1), 
                                                               _mm512_set1_pd(1.737331760720576030932e-8)};
	                   __ATTR_ALIGN__(64) static const __m512d d[9]  =  {_mm512_set1_pd(8.258160008564488034698e-2), 
	                                                        _mm512_set1_pd(4.344836335509282083360e+0), 
                                                                _mm512_set1_pd(4.662179610356861756812e+1), 
                                                                _mm512_set1_pd(1.775728186717289799677e+2), 
                                                                _mm512_set1_pd(2.953136335677908517423e+2), 
                                                                _mm512_set1_pd(2.342573504717625153053e+2), 
                                                                _mm512_set1_pd(9.021658450529372642314e+1), 
                                                                _mm512_set1_pd(1.587964570758947927903e+1), 
                                                                _mm512_set1_pd(1.000000000000000000000e+0)};                                  
	                   __ATTR_ALIGN__(64) static const __m512d e[10] = {_mm512_set1_pd(1.1669552669734461083368e+2),
	                                                       _mm512_set1_pd(2.1500672908092918123209e+3),
	                                                       _mm512_set1_pd(1.5924175980637303639884e+4),
	                                                       _mm512_set1_pd(8.9904972007457256553251e+4),
	                                                       _mm512_set1_pd(1.5026059476436982420737e+5),
	                                                       _mm512_set1_pd(-1.4815102102575750838086e+5),
	                                                       _mm512_set1_pd(5.0196785185439843791020e+0)};
	                   __ATTR_ALIGN__(64) static const __m512d f[10] = {_mm512_set1_pd(3.9147856245556345627078e+4),
	                                                       _mm512_set1_pd(2.5989762083608489777411e+5),
	                                                       _mm512_set1_pd(5.5903756210022864003380e+5),
	                                                       _mm512_set1_pd(5.4616842050691155735758e+5),
	                                                       _mm512_set1_pd(2.7858134710520842139357e+5),
	                                                       _mm512_set1_pd(7.9231787945279043698718e+4),
	                                                       _mm512_set1_ps(1.2842808586627297365998e+4),
	                                                       _mm512_set1_ps(1.1635769915320848035459e+3),
	                                                       _mm512_set1_ps(5.4199632588522559414924e+1),
	                                                       _mm512_set1_ps(1.0)};
	                  __ATTR_ALIGN__(64) static const __m512d plg[4] = {_mm512_set1_pd(2.4562334077563243311e+01),
	                                                       _mm512_set1_pd(2.3642701335621505212e+02), 
                                                               _mm512_set1_pd(-5.4989956895857911039e+02),
                                                               _mm512_set1_pd(3.5687548468071500413e+02)};
	                  __ATTR_ALIGN__(64) static const __m512d qlg[4] = {_mm512_set1_pd(-3.5553900764052419184e+01),
	                                                       _mm512_set1_pd(1.9400230218539473193e+02), 
                                                               _mm512_set1_pd(-3.3442903192607538956e+02),
                                                               _mm512_set1_pd(1.7843774234035750207e+02)};
	                   __ATTR_ALIGN__(64) static const __m512d p[10]  = {_mm512_set1_pd(-1.2963702602474830028590e+1),
	                                                        _mm512_set1_pd(-1.2831220659262000678155e+3),
                                                                _mm512_set1_pd(-1.4287072500197005777376e+4),
                                                                _mm512_set1_pd(-1.4299841572091610380064e+6), 
                                                                _mm512_set1_pd(-3.1398660864247265862050e+5),
                                                                _mm512_set1_pd(-3.5377809694431133484800e+8),
                                                                _mm512_set1_pd(3.1984354235237738511048e+8),
                                                                _mm512_set1_pd(-2.5301823984599019348858e+10),
                                                                _mm512_set1_pd(1.2177698136199594677580e+10),
                                                                _mm512_set1_pd(-2.0829040666802497120940e+11)};
                           __ATTR_ALIGN__(64) static const __m512d q[10]  = {_mm512_set1_pd(7.6886718750000000000000e+1),
                                                                _mm512_set1_pd(-5.5648470543369082846819e+3),
                                                                _mm512_set1_pd(1.9418469440759880361415e+5),
                                                                _mm512_set1_pd(-4.2648434812177161405483e+6), 
                                                                _mm512_set1_pd(6.4698830956576428587653e+7),
                                                                _mm512_set1_pd(-7.0108568774215954065376e+8), 
                                                                _mm512_set1_pd(5.4229617984472955011862e+9),
                                                                _mm512_set1_pd(-2.8986272696554495342658e+10,
                                                                _mm512_set1_pd(9.8900934262481749439886e+10),
                                                                _mm512_set1_pd(-8.9673749185755048616855e+10)};
	                   __ATTR_ALIGN__(64) static const __m512d p1[10] = {_mm512_set1_pd(1.647721172463463140042e+0),
	                                                        _mm512_set1_pd(-1.860092121726437582253e+1),
                                                                _mm512_set1_pd(-1.000641913989284829961e+1),
                                                                _mm512_set1_pd(-2.105740799548040450394e+1), 
                                                                _mm512_set1_pd(-9.134835699998742552432e-1),
                                                                _mm512_set1_pd(-3.323612579343962284333e+1),
                                                                _mm512_set1_pd(2.495487730402059440626e+1), 
                                                                _mm512_set1_pd(2.652575818452799819855e+1), 
                                                                _mm512_set1_pd(-1.845086232391278674524e+0), 
                                                                _mm512_set1_pd(9.999933106160568739091e-1)};
                           __ATTR_ALIGN__(64) static const __m512d q1[9]  = {_mm512_set1_pd(9.792403599217290296840e+1), 
                                                                _mm512_set1_pd(6.403800405352415551324e+1),
                                                                _mm512_set1_pd(5.994932325667407355255e+1), 
                                                                _mm512_set1_pd(2.538819315630708031713e+2),
                                                                _mm512_set1_pd(4.429413178337928401161e+1), 
                                                                _mm512_set1_pd(1.192832423968601006985e+3), 
                                                                _mm512_set1_pd(1.991004470817742470726e+2),
                                                                _mm512_set1_pd(-1.093556195391091143924e+1), 
                                                                _mm512_set1_pd(1.001533852045342697818e+0)};
	                   __ATTR_ALIGN__(64) static const __m512d p2[10]  = {_mm512_set1_pd(1.75338801265465972390e+2),
	                                                         _mm512_set1_pd(-2.23127670777632409550e+2), 
                                                                 _mm512_set1_pd(-1.81949664929868906455e+1),
                                                                 _mm512_set1_pd(-2.79798528624305389340e+1), 
                                                                 _mm512_set1_pd(-7.63147701620253630855e+0),
                                                                 _mm512_set1_pd(-1.52856623636929636839e+1), 
                                                                 _mm512_set1_pd(-7.06810977895029358836e+0),
                                                                 _mm512_set1_pd(-5.00006640413131002475e+0),
                                                                 _mm512_set1_pd(-3.00000000320981265753e+0), 
                                                                 _mm512_set1_pd(1.00000000000000485503e+0)};
                            __ATTR_ALIGN__(64) static const __m512d q2[9]  = {_mm512_set1_pd(3.97845977167414720840e+4), 
                                                                 _mm512_set1_pd(3.97277109100414518365e+0), 
                                                                 _mm512_set1_pd(1.37790390235747998793e+2), 
                                                                 _mm512_set1_pd(1.17179220502086455287e+2), 
                                                                 _mm512_set1_pd(7.04831847180424675988e+1),
                                                                 _mm512_set1_pd(-1.20187763547154743238e+1), 
                                                                 _mm512_set1_pd(-7.99243595776339741065e+0),
                                                                 _mm512_set1_pd(-2.99999894040324959612e+0),
                                                                 _mm512_set1_pd(1.99999999999048104167e+0)};   
                            __ATTR_ALIGN__(64) static const __m512d r[10]  = {_mm512_set1_pd(2.645677793077147237806e+0),
                                                                 _mm512_set1_pd(-2.378372882815725244124e+0), 
                                                                 _mm512_set1_pd(-2.421106956980653511550e+1), 
                                                                 _mm512_set1_pd(1.052976392459015155422e+1), 
                                                                 _mm512_set1_pd(1.945603779539281810439e+1),
                                                                 _mm512_set1_pd(-3.015761863840593359165e+1), 
                                                                 _mm512_set1_pd(1.120011024227297451523e+1),
                                                                 _mm512_set1_pd(-3.988850730390541057912e+0), 
                                                                 _mm512_set1_pd(9.565134591978630774217e+0), 
                                                                 _mm512_set1_pd(9.981193787537396413219e-1)};
                           __ATTR_ALIGN__(64) static const __m512d s[9]    = {_mm512_set1_pd(1.598517957704779356479e-4), 
                                                                 _mm512_set1_pd(4.644185932583286942650e+0), 
                                                                 _mm512_set1_pd(3.697412299772985940785e+2),
                                                                 _mm512_set1_pd(-8.791401054875438925029e+0), 
                                                                 _mm512_set1_pd(7.608194509086645763123e+2), 
                                                                 _mm512_set1_pd(2.852397548119248700147e+1), 
                                                                 _mm512_set1_pd(4.731097187816050252967e+2),
                                                                 _mm512_set1_pd(-2.369210235636181001661e+2), 
                                                                 _mm512_set1_pd(1.249884822712447891440e+0)};  
                           __ATTR_ALIGN__(64) __m512d q[10];
	                   __ATTR_ALIGN__(64) __m512d qlq[10];
	                   __ATTR_ALIGN__(64) __m512d qx[10];
	                   __ATTR_ALIGN__(64) __m512d px[10];                                       
                          const __m512d zero = _mm512_set1_pd(0.0e+0);
                          const __m512d p037 = _mm512_set1_pd(0.037e+0);
                          const __m512d half = _mm512_set1_pd(0.5);
                          const __m512d one  = _mm512_set1_pd(1.0e+0);
                          const __m512d two  = _mm512_set1_pd(2.0e+0);
                          const __m512d three= _mm512_set1_pd(3.0e+0);
                          const __m512d four = _mm512_set1_pd(4.0e+0);
                          const __m512d six  = _mm512_set1_pd(6.0e+0);
                          const __m512d twelve = _mm512_set1_pd(2.0e+0);
                          const __m512d two4 = _mm512_set1_pd(24.0e+0);
                          const __m512d fourty = _mm512_set1_pd(40.0e+0);
                          const __m512d exp40 = _mm512_set1_pd(2.3538526683701998541e+17);
                          const __m512d x01  = _mm512_set1_pd(381.5e+0);
                          const __m512d x11  = _mm512_set1_pd(1024.0e+0);
                          const __m512d x02  = _mm512_set1_pd(-5.1182968633365538008e-5);
                          const __m512d x0   = _mm512_set1_pd(3.7250741078136663466e-1);   
                          const __m512d xinf = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax = _mm512_set1_pd(716.351e+0);
                          const __m512d xbig = _mm512_set1_pd(701.84e+0);                        
	                  __m512d ei,frac,result,sump,sumq;
	                  __m512d t,w,x;
	                  __m512d mx0,y,ysq;
	                 
	                  x = arg;
	                  if(__m512_cmp_pd_mask(x,zero,_CMP_EQ_OQ)) {
	                  
	                      ei = negate_zmm8r8(xinf);
	                      if(jint == 2) {
	                      
	                          ei = negate_zmm8r8(ei);
	                      }
	                     /*
	                        !
                                 !  Calculate EI for negative argument or for E1.
                                !   
	                     */
	                 }
	                 else if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ) || 
	                         jint == 2) {
	                         
	                         y = _mm512_abs_pd(x);
	                         if(_mm512_cmp_pd_mask(y,one,_CMP_LT)) {
	                         
	                            sump = _mm512_fmadd_pd(a[6],y,a[0]);
	                            sumq = _mm512_add_pd(y,b[0]);
	                            
	                            sump = _mm512_fmadd_pd(sump,y,a[1]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[1]);
	                            sump = _mm512_fmadd_pd(sump,y,a[2]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[2]);
	                            sump = _mm512_fmadd_pd(sump,y,a[3]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[3]);
	                            sump = _mm512_fmadd_pd(sump,y,a[4]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[4]);
	                            sump = _mm512_fmadd_pd(sump,y,a[5]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[5]);
	                            
	                            ei   = _mm512_sub_pd(xlog(y),
	                                             _mm512_div_pd(sump,sumq));
	                            if(jint==3) {
	                               ei = _mm512_mul_pd(ei,xexp(y);
	                            }
	                            
	                         }  
	                         else if(_mm512_cmp_pd_mask(y,four,_CMP_LE_OQ)) {
	                        
	                                w    = _mm512_div_pd(one,y);
	                                sump = c[0];
	                                sumq = d[0];
	                                
	                                sump = _mm512_fmadd_pd(sump,w,c[1]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[1]);
	                                sump = _mm512_fmadd_pd(sump,w,c[2]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[2]);
	                                sump = _mm512_fmadd_pd(sump,w,c[3]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[3]);  
	                                sump = _mm512_fmadd_pd(sump,w,c[4]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[4]);
	                                sump = _mm512_fmadd_pd(sump,w,c[5]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[5]);
	                                sump = _mm512_fmadd_pd(sump,w,c[6]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[6]);
	                                sump = _mm512_fmadd_pd(sump,w,c[7]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[7]);
	                                sump = _mm512_fmadd_pd(sump,w,c[8]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[8]);
	                                
	                                ei   = _mm512_div_pd(negate_zmm8r8(sump),sumq);
	                                
	                                if(jint!=3) {
	                                    ei = _mm512_mul_pd(ei,
	                                                    xexp(negate_zmm8r8(y)));
	                                }
	                        } 
	                        else {
	                                 
	                                if(_mm512_cmp_pd_mask(xbig,y,_CMP_LT_OQ) &&
	                                   jint<3) {
	                                     ei = zero;  
	                                }
                                        else {
                                           
                                             w    = _mm512_div_pd(one,y);
                                             sump = e[0];
                                             sumq = f[0];
                                             
                                             sump = _mm512_fmadd_pd(sump,w,e[1]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[1]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[2]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[2]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[3]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[3]);  
	                                     sump = _mm512_fmadd_pd(sump,w,e[4]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[4]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[5]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[5]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[6]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[6]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[7]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[7]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[8]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[8]); 
	                                     sump = _mm512_fmadd_pd(sump,w,e[9]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[9]); 
	                                     
	                                     const __m512 x0 = _mm512_div_pd(sump,sumq);
	                                     const __m512 x1 = _mm512_sub_pd(one,w);
	                                     ei    = _mm512_mul_pd(negate_zmm8r8(w),
	                                                       _mm512_mul_pd(x0,x1));
	                                                       
	                                     if(jint!=3) {
	                                         ei = _mm512_mul_pd(ei,
	                                                         negate_zmm8r8(y));
	                                     }
	                                }	                                
	                                
	                        }    
	                        
	                        if(jint==2) {
	                           ei = negate_zmm8r8(ei);
	                        } 
	                        
	                        /*
	                           !
                                   !  To improve conditioning, rational approximations are expressed
                                   !  in terms of Chebyshev polynomials for 0 <= X < 6, and in
                                   !  continued fraction form for larger X.
                                   !
	                        */
	                        
	                }
	                else if(_mm512_cmp_pd_mask(x,six,_CMP_LT_OQ)) {
	                       
	                         t     = _mm512_add_pd(x,x);
	                         t     = _mm512_sub_pd(_mm512_div_pd(t,three),two);
	                         px[0] = zero;
	                         qx[0] = zero;
	                         px[1] = p[0];
	                         qx[1] = q[0];
	                         
	                         px[2] = _mm512_fmsub_pd(t,px[1],_mm512_add_pd(px[0],p[1]));
	                         qx[2] = _mm512_fmsub_pd(t,qx[1],_mm512_add_pd(qx[0],q[1]));
	                         px[3] = _mm512_fmsub_pd(t,px[2],_mm512_add_pd(px[1],p[2]));
	                         qx[3] = _mm512_fmsub_pd(t,qx[2],_mm512_add_pd(qx[1],q[2]));
	                         px[4] = _mm512_fmsub_pd(t,px[3],_mm512_add_pd(px[2],p[3]));
	                         qx[4] = _mm512_fmsub_pd(t,qx[3],_mm512_add_pd(qx[2],q[3]));
	                         px[5] = _mm512_fmsub_pd(t,px[4],_mm512_add_pd(px[3],p[4]));
	                         qx[5] = _mm512_fmsub_pd(t,qx[4],_mm512_add_pd(qx[3],q[4]));
	                         px[6] = _mm512_fmsub_pd(t,px[5],_mm512_add_pd(px[4],p[5]));
	                         qx[6] = _mm512_fmsub_pd(t,qx[5],_mm512_add_pd(qx[4],q[5]));
	                         px[7] = _mm512_fmsub_pd(t,px[6],_mm512_add_pd(px[5],p[6]));
	                         qx[7] = _mm512_fmsub_pd(t,qx[6],_mm512_add_pd(qx[5],q[6]));
	                         px[8] = _mm512_fmsub_pd(t,px[7],_mm512_add_pd(px[6],p[7]));
	                         qx[8] = _mm512_fmsub_pd(t,qx[7],_mm512_add_pd(qx[6],q[7]));
	                         __m512 tmp = _mm512_mul_pd(half,t);
	                         sump  = _mm512_fmsub_pd(tmp,px[9],_mm512_add_pd(px[8],p[9]));
	                         sumq  = _mm512_fmsub_pd(tmp,qx[9],_mm512_add_pd(qx[8],q[9]));
	                         frac  = _mm512_div_pd(sump,sumq);
	                         tmp    = _mm512_div_sub_pd(x,_mm512_div_pd(x01,x11));
	                         xmx0  = _mm512_sub_pd(tmp,x02);

                                 if(_mm512_cmp_pd_mask(p037,
                                                   _mm512_abs_pd(xmx0),_CMP_LE_OQ)) {
                                     
                                     __m512 tmp = _mm512_div_pd(x,x0);
                                     ei         = _mm512_fmadd_pd(frac,xmx0,xlog(tmp));
                                     if(jint==3) {
                                        ei = _mm512_mul_pd(xexp(negate_zmm8r8(x),ei));
                                     }                 
                                 }
                                 else {
                                    
                                      //Special approximation to ln(X/X0) for X close to X0. 
                                      
                                      y    = _mm512_div_pd(xmx0,_mm512_add_pd(x,x0));
                                      ysq  = _mm512_mul_pd(y,y);
                                      sump = plg[0]; 
                                      sumq = _mm512_add_pd(ysq,qlg[0]);
                                      
                                      sump = _mm512_fmadd_pd(sump,ysq,plg[1]);
                                      sumq = _mm512_fmadd_pd(sumq,ysq,qlg[1]);
                                      sump = _mm512_fmadd_pd(sump,ysq,plg[2]);
                                      sumq = _mm512_fmadd_pd(sumq,ysq,qlg[2]);
                                      sump = _mm512_fmadd_pd(sump,ysq,plg[3]);
                                      sumq = _mm512_fmadd_pd(sumq,ysq,qlg[3]);
                                      
                                      __m512 tmp = _mm512_fmadd_pd(sumq,
                                                               _mm512_add_pd(x,x0),frac);
                                      ei   = _mm512_mul_pd(_mm512_div_pd(sump,tmp),xmx0);
                                      
                                      if(jint==3) {
                                         ei = _mm512_mul_pd(xexp(negate_zmm8r8(x),ei));
                                      }
                                 }
	                }
	                else if(_mm512_cmp_pd_mask(x,twelve,_CMP_LT_OQ)) {
	                         
	                         frac = zero;
	                         
	                         frac = _mm512_div_pd(s[0],
	                                          _mm512_add_pd(r[0],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[1],
	                                          _mm512_add_pd(r[1],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[2],
	                                          _mm512_add_pd(r[2],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[3],
	                                          _mm512_add_pd(r[3],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[4],
	                                          _mm512_add_pd(r[4],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[5],
	                                          _mm512_add_pd(r[5],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[6],
	                                          _mm512_add_pd(r[6],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[7],
	                                          _mm512_add_pd(r[7],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[8],
	                                          _mm512_add_pd(r[8],
	                                                    _mm512_add_pd(frac)));
	                                                    
	                         ei   = _mm512_div_pd(_mm512_add_pd(r[9],frac),x);
	                         if(jint!=3) {
	                            ei = _mm512_mul_pd(ei,xexp(x));
	                         }
	                         
	                }
	                else if(_mm512_cmp_pd_mask(x,two4,_CMP_LE_OQ)) {
	                         
	                         frac = zero;
	                         
	                         frac = _mm512_div_pd(q1[0],
	                                          _mm512_add_pd(p1[0],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[1],
	                                          _mm512_add_pd(p1[1],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[2],
	                                          _mm512_add_pd(p1[2],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[3],
	                                          _mm512_add_pd(p1[3],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[4],
	                                          _mm512_add_pd(p1[4],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[5],
	                                          _mm512_add_pd(p1[5],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[6],
	                                          _mm512_add_pd(p1[6],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[7],
	                                          _mm512_add_pd(p1[7],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[8],
	                                          _mm512_add_pd(p1[8],
	                                                    _mm512_add_pd(frac)));
	                                                    
	                         ei   = _mm512_div_pd(_mm512_add_pd(p1[9],frac),x);
	                         if(jint!=3) {
	                            ei = _mm512_mul_pd(ei,xexp(x));
	                         }
	                         
	                }
	                else {
	                   
	                       if(_mm512_cmp_pd_mask(xmax,x,_CMP_LE_OQ) &&
	                          jint==3) {
	                          
	                           ei = xinf;  
	                       }
	                       else {
	                          
	                           y    = _mm512_div_pd(one,x);
	                           frac = zero;
	                           
	                           frac = _mm512_div_pd(q2[0],
	                                          _mm512_add_pd(p2[0],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[1],
	                                          _mm512_add_pd(p2[1],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[2],
	                                          _mm512_add_pd(p2[2],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[3],
	                                          _mm512_add_pd(p2[3],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[4],
	                                          _mm512_add_pd(p2[4],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[5],
	                                          _mm512_add_pd(p2[5],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[6],
	                                          _mm512_add_pd(p2[6],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[7],
	                                          _mm512_add_pd(p2[7],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[8],
	                                          _mm512_add_pd(p2[8],
	                                                    _mm512_add_pd(frac)));
	                                                    
	                           frac = _mm512_add_pd(p2[9],frac);
	                           ei   = _mm512_fmadd_pd(frac,
	                                              _mm512_mul_pd(y,y),y);
	                           if(jint!=3) {
	                              
	                              if(_mm512_cmp_pd_mask(x,
	                                       _mm512_sub_pd(xmax,two4),_CMP_LE_OQ)) {
	                                  
	                                  ei = _mm512_mul_pd(ei,xexp(x));    
	                                  // Calculation reformulated to avoid premature overflow.    
	                              }
	                              else {
	                                 
	                                  __m512 tmp = _mm512_sub_pd(x,fourty);
	                                  ei  = _mm512_mul_pd(_mm512_mul_pd(ei,xexp(tmp),exp40));
	                              }
	                           }
	                       }
	                }
	                 
	                 result = ei;
	                 return (result);
	                                
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calcei_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                   const int32_t jint) {
	                                 
	                                 
	                   __ATTR_ALIGN__(64) static const __m512d a[7]  = {_mm512_set1_pd(1.1669552669734461083368e+2),
	                                                       _mm512_set1_pd(2.1500672908092918123209e+3), 
                                                               _mm512_set1_pd(1.5924175980637303639884e+4), 
                                                               _mm512_set1_pd(8.9904972007457256553251e+4), 
                                                               _mm512_set1_pd(1.5026059476436982420737e+5),
                                                               _mm512_set1_pd(-1.4815102102575750838086e+5), 
                                                               _mm512_set1_pd(5.0196785185439843791020e+0)};
                           __ATTR_ALIGN__(64) static const __m512d b[6]  = {_mm512_set1_pd(4.0205465640027706061433e+1), 
                                                               _mm512_set1_pd(7.5043163907103936624165e+2),
                                                               _mm512_set1_pd(8.1258035174768735759855e+3), 
                                                               _mm512_set1_pd(5.2440529172056355429883e+4), 
                                                               _mm512_set1_pd(1.8434070063353677359298e+5), 
                                                               _mm512_set1_pd(2.5666493484897117319268e+5)};
                           __ATTR_ALIGN__(64) static const __m512d c[9]  = {_mm512_set1_pd(3.828573121022477169108e-1), 
                                                               _mm512_set1_pd(1.107326627786831743809e+1), 
                                                               _mm512_set1_pd(7.246689782858597021199e+1), 
                                                               _mm512_set1_pd(1.700632978311516129328e+2), 
                                                               _mm512_set1_pd(1.698106763764238382705e+2), 
                                                               _mm512_set1_pd(7.633628843705946890896e+1), 
                                                               _mm512_set1_pd(1.487967702840464066613e+1), 
                                                               _mm512_set1_pd(9.999989642347613068437e-1), 
                                                               _mm512_set1_pd(1.737331760720576030932e-8)};
	                   __ATTR_ALIGN__(64) static const __m512d d[9]  =  {_mm512_set1_pd(8.258160008564488034698e-2), 
	                                                        _mm512_set1_pd(4.344836335509282083360e+0), 
                                                                _mm512_set1_pd(4.662179610356861756812e+1), 
                                                                _mm512_set1_pd(1.775728186717289799677e+2), 
                                                                _mm512_set1_pd(2.953136335677908517423e+2), 
                                                                _mm512_set1_pd(2.342573504717625153053e+2), 
                                                                _mm512_set1_pd(9.021658450529372642314e+1), 
                                                                _mm512_set1_pd(1.587964570758947927903e+1), 
                                                                _mm512_set1_pd(1.000000000000000000000e+0)};                                  
	                   __ATTR_ALIGN__(64) static const __m512d e[10] = {_mm512_set1_pd(1.1669552669734461083368e+2),
	                                                       _mm512_set1_pd(2.1500672908092918123209e+3),
	                                                       _mm512_set1_pd(1.5924175980637303639884e+4),
	                                                       _mm512_set1_pd(8.9904972007457256553251e+4),
	                                                       _mm512_set1_pd(1.5026059476436982420737e+5),
	                                                       _mm512_set1_pd(-1.4815102102575750838086e+5),
	                                                       _mm512_set1_pd(5.0196785185439843791020e+0)};
	                   __ATTR_ALIGN__(64) static const __m512d f[10] = {_mm512_set1_pd(3.9147856245556345627078e+4),
	                                                       _mm512_set1_pd(2.5989762083608489777411e+5),
	                                                       _mm512_set1_pd(5.5903756210022864003380e+5),
	                                                       _mm512_set1_pd(5.4616842050691155735758e+5),
	                                                       _mm512_set1_pd(2.7858134710520842139357e+5),
	                                                       _mm512_set1_pd(7.9231787945279043698718e+4),
	                                                       _mm512_set1_ps(1.2842808586627297365998e+4),
	                                                       _mm512_set1_ps(1.1635769915320848035459e+3),
	                                                       _mm512_set1_ps(5.4199632588522559414924e+1),
	                                                       _mm512_set1_ps(1.0)};
	                  __ATTR_ALIGN__(64) static const __m512d plg[4] = {_mm512_set1_pd(2.4562334077563243311e+01),
	                                                       _mm512_set1_pd(2.3642701335621505212e+02), 
                                                               _mm512_set1_pd(-5.4989956895857911039e+02),
                                                               _mm512_set1_pd(3.5687548468071500413e+02)};
	                  __ATTR_ALIGN__(64) static const __m512d qlg[4] = {_mm512_set1_pd(-3.5553900764052419184e+01),
	                                                       _mm512_set1_pd(1.9400230218539473193e+02), 
                                                               _mm512_set1_pd(-3.3442903192607538956e+02),
                                                               _mm512_set1_pd(1.7843774234035750207e+02)};
	                   __ATTR_ALIGN__(64) static const __m512d p[10]  = {_mm512_set1_pd(-1.2963702602474830028590e+1),
	                                                        _mm512_set1_pd(-1.2831220659262000678155e+3),
                                                                _mm512_set1_pd(-1.4287072500197005777376e+4),
                                                                _mm512_set1_pd(-1.4299841572091610380064e+6), 
                                                                _mm512_set1_pd(-3.1398660864247265862050e+5),
                                                                _mm512_set1_pd(-3.5377809694431133484800e+8),
                                                                _mm512_set1_pd(3.1984354235237738511048e+8),
                                                                _mm512_set1_pd(-2.5301823984599019348858e+10),
                                                                _mm512_set1_pd(1.2177698136199594677580e+10),
                                                                _mm512_set1_pd(-2.0829040666802497120940e+11)};
                           __ATTR_ALIGN__(64) static const __m512d q[10]  = {_mm512_set1_pd(7.6886718750000000000000e+1),
                                                                _mm512_set1_pd(-5.5648470543369082846819e+3),
                                                                _mm512_set1_pd(1.9418469440759880361415e+5),
                                                                _mm512_set1_pd(-4.2648434812177161405483e+6), 
                                                                _mm512_set1_pd(6.4698830956576428587653e+7),
                                                                _mm512_set1_pd(-7.0108568774215954065376e+8), 
                                                                _mm512_set1_pd(5.4229617984472955011862e+9),
                                                                _mm512_set1_pd(-2.8986272696554495342658e+10,
                                                                _mm512_set1_pd(9.8900934262481749439886e+10),
                                                                _mm512_set1_pd(-8.9673749185755048616855e+10)};
	                   __ATTR_ALIGN__(64) static const __m512d p1[10] = {_mm512_set1_pd(1.647721172463463140042e+0),
	                                                        _mm512_set1_pd(-1.860092121726437582253e+1),
                                                                _mm512_set1_pd(-1.000641913989284829961e+1),
                                                                _mm512_set1_pd(-2.105740799548040450394e+1), 
                                                                _mm512_set1_pd(-9.134835699998742552432e-1),
                                                                _mm512_set1_pd(-3.323612579343962284333e+1),
                                                                _mm512_set1_pd(2.495487730402059440626e+1), 
                                                                _mm512_set1_pd(2.652575818452799819855e+1), 
                                                                _mm512_set1_pd(-1.845086232391278674524e+0), 
                                                                _mm512_set1_pd(9.999933106160568739091e-1)};
                           __ATTR_ALIGN__(64) static const __m512d q1[9]  = {_mm512_set1_pd(9.792403599217290296840e+1), 
                                                                _mm512_set1_pd(6.403800405352415551324e+1),
                                                                _mm512_set1_pd(5.994932325667407355255e+1), 
                                                                _mm512_set1_pd(2.538819315630708031713e+2),
                                                                _mm512_set1_pd(4.429413178337928401161e+1), 
                                                                _mm512_set1_pd(1.192832423968601006985e+3), 
                                                                _mm512_set1_pd(1.991004470817742470726e+2),
                                                                _mm512_set1_pd(-1.093556195391091143924e+1), 
                                                                _mm512_set1_pd(1.001533852045342697818e+0)};
	                   __ATTR_ALIGN__(64) static const __m512d p2[10]  = {_mm512_set1_pd(1.75338801265465972390e+2),
	                                                         _mm512_set1_pd(-2.23127670777632409550e+2), 
                                                                 _mm512_set1_pd(-1.81949664929868906455e+1),
                                                                 _mm512_set1_pd(-2.79798528624305389340e+1), 
                                                                 _mm512_set1_pd(-7.63147701620253630855e+0),
                                                                 _mm512_set1_pd(-1.52856623636929636839e+1), 
                                                                 _mm512_set1_pd(-7.06810977895029358836e+0),
                                                                 _mm512_set1_pd(-5.00006640413131002475e+0),
                                                                 _mm512_set1_pd(-3.00000000320981265753e+0), 
                                                                 _mm512_set1_pd(1.00000000000000485503e+0)};
                            __ATTR_ALIGN__(64) static const __m512d q2[9]  = {_mm512_set1_pd(3.97845977167414720840e+4), 
                                                                 _mm512_set1_pd(3.97277109100414518365e+0), 
                                                                 _mm512_set1_pd(1.37790390235747998793e+2), 
                                                                 _mm512_set1_pd(1.17179220502086455287e+2), 
                                                                 _mm512_set1_pd(7.04831847180424675988e+1),
                                                                 _mm512_set1_pd(-1.20187763547154743238e+1), 
                                                                 _mm512_set1_pd(-7.99243595776339741065e+0),
                                                                 _mm512_set1_pd(-2.99999894040324959612e+0),
                                                                 _mm512_set1_pd(1.99999999999048104167e+0)};   
                            __ATTR_ALIGN__(64) static const __m512d r[10]  = {_mm512_set1_pd(2.645677793077147237806e+0),
                                                                 _mm512_set1_pd(-2.378372882815725244124e+0), 
                                                                 _mm512_set1_pd(-2.421106956980653511550e+1), 
                                                                 _mm512_set1_pd(1.052976392459015155422e+1), 
                                                                 _mm512_set1_pd(1.945603779539281810439e+1),
                                                                 _mm512_set1_pd(-3.015761863840593359165e+1), 
                                                                 _mm512_set1_pd(1.120011024227297451523e+1),
                                                                 _mm512_set1_pd(-3.988850730390541057912e+0), 
                                                                 _mm512_set1_pd(9.565134591978630774217e+0), 
                                                                 _mm512_set1_pd(9.981193787537396413219e-1)};
                           __ATTR_ALIGN__(64) static const __m512d s[9]    = {_mm512_set1_pd(1.598517957704779356479e-4), 
                                                                 _mm512_set1_pd(4.644185932583286942650e+0), 
                                                                 _mm512_set1_pd(3.697412299772985940785e+2),
                                                                 _mm512_set1_pd(-8.791401054875438925029e+0), 
                                                                 _mm512_set1_pd(7.608194509086645763123e+2), 
                                                                 _mm512_set1_pd(2.852397548119248700147e+1), 
                                                                 _mm512_set1_pd(4.731097187816050252967e+2),
                                                                 _mm512_set1_pd(-2.369210235636181001661e+2), 
                                                                 _mm512_set1_pd(1.249884822712447891440e+0)};  
                           __ATTR_ALIGN__(64) __m512d q[10];
	                   __ATTR_ALIGN__(64) __m512d qlq[10];
	                   __ATTR_ALIGN__(64) __m512d qx[10];
	                   __ATTR_ALIGN__(64) __m512d px[10];                                       
                          const __m512d zero = _mm512_set1_pd(0.0e+0);
                          const __m512d p037 = _mm512_set1_pd(0.037e+0);
                          const __m512d half = _mm512_set1_pd(0.5);
                          const __m512d one  = _mm512_set1_pd(1.0e+0);
                          const __m512d two  = _mm512_set1_pd(2.0e+0);
                          const __m512d three= _mm512_set1_pd(3.0e+0);
                          const __m512d four = _mm512_set1_pd(4.0e+0);
                          const __m512d six  = _mm512_set1_pd(6.0e+0);
                          const __m512d twelve = _mm512_set1_pd(2.0e+0);
                          const __m512d two4 = _mm512_set1_pd(24.0e+0);
                          const __m512d fourty = _mm512_set1_pd(40.0e+0);
                          const __m512d exp40 = _mm512_set1_pd(2.3538526683701998541e+17);
                          const __m512d x01  = _mm512_set1_pd(381.5e+0);
                          const __m512d x11  = _mm512_set1_pd(1024.0e+0);
                          const __m512d x02  = _mm512_set1_pd(-5.1182968633365538008e-5);
                          const __m512d x0   = _mm512_set1_pd(3.7250741078136663466e-1);   
                          const __m512d xinf = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax = _mm512_set1_pd(716.351e+0);
                          const __m512d xbig = _mm512_set1_pd(701.84e+0);                        
	                  __m512d arg,ei,frac,result,sump,sumq;
	                  __m512d t,w,x;
	                  __m512d mx0,y,ysq;
	                 
	                  arg = _mm512_load_pd(&parg[0]);
	                  x = arg;
	                  if(__m512_cmp_pd_mask(x,zero,_CMP_EQ_OQ)) {
	                  
	                      ei = negate_zmm8r8(xinf);
	                      if(jint == 2) {
	                      
	                          ei = negate_zmm8r8(ei);
	                      }
	                     /*
	                        !
                                 !  Calculate EI for negative argument or for E1.
                                !   
	                     */
	                 }
	                 else if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ) || 
	                         jint == 2) {
	                         
	                         y = _mm512_abs_pd(x);
	                         if(_mm512_cmp_pd_mask(y,one,_CMP_LT)) {
	                         
	                            sump = _mm512_fmadd_pd(a[6],y,a[0]);
	                            sumq = _mm512_add_pd(y,b[0]);
	                            
	                            sump = _mm512_fmadd_pd(sump,y,a[1]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[1]);
	                            sump = _mm512_fmadd_pd(sump,y,a[2]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[2]);
	                            sump = _mm512_fmadd_pd(sump,y,a[3]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[3]);
	                            sump = _mm512_fmadd_pd(sump,y,a[4]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[4]);
	                            sump = _mm512_fmadd_pd(sump,y,a[5]);
	                            sumq = _mm512_fmadd_pd(sumq,y,b[5]);
	                            
	                            ei   = _mm512_sub_pd(xlog(y),
	                                             _mm512_div_pd(sump,sumq));
	                            if(jint==3) {
	                               ei = _mm512_mul_pd(ei,xexp(y);
	                            }
	                            
	                         }  
	                         else if(_mm512_cmp_pd_mask(y,four,_CMP_LE_OQ)) {
	                        
	                                w    = _mm512_div_pd(one,y);
	                                sump = c[0];
	                                sumq = d[0];
	                                
	                                sump = _mm512_fmadd_pd(sump,w,c[1]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[1]);
	                                sump = _mm512_fmadd_pd(sump,w,c[2]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[2]);
	                                sump = _mm512_fmadd_pd(sump,w,c[3]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[3]);  
	                                sump = _mm512_fmadd_pd(sump,w,c[4]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[4]);
	                                sump = _mm512_fmadd_pd(sump,w,c[5]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[5]);
	                                sump = _mm512_fmadd_pd(sump,w,c[6]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[6]);
	                                sump = _mm512_fmadd_pd(sump,w,c[7]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[7]);
	                                sump = _mm512_fmadd_pd(sump,w,c[8]);
	                                sumq = _mm512_fmadd_pd(sumq,w,d[8]);
	                                
	                                ei   = _mm512_div_pd(negate_zmm8r8(sump),sumq);
	                                
	                                if(jint!=3) {
	                                    ei = _mm512_mul_pd(ei,
	                                                    xexp(negate_zmm8r8(y)));
	                                }
	                        } 
	                        else {
	                                 
	                                if(_mm512_cmp_pd_mask(xbig,y,_CMP_LT_OQ) &&
	                                   jint<3) {
	                                     ei = zero;  
	                                }
                                        else {
                                           
                                             w    = _mm512_div_pd(one,y);
                                             sump = e[0];
                                             sumq = f[0];
                                             
                                             sump = _mm512_fmadd_pd(sump,w,e[1]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[1]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[2]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[2]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[3]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[3]);  
	                                     sump = _mm512_fmadd_pd(sump,w,e[4]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[4]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[5]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[5]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[6]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[6]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[7]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[7]);
	                                     sump = _mm512_fmadd_pd(sump,w,e[8]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[8]); 
	                                     sump = _mm512_fmadd_pd(sump,w,e[9]);
	                                     sumq = _mm512_fmadd_pd(sumq,w,f[9]); 
	                                     
	                                     const __m512 x0 = _mm512_div_pd(sump,sumq);
	                                     const __m512 x1 = _mm512_sub_pd(one,w);
	                                     ei    = _mm512_mul_pd(negate_zmm8r8(w),
	                                                       _mm512_mul_pd(x0,x1));
	                                                       
	                                     if(jint!=3) {
	                                         ei = _mm512_mul_pd(ei,
	                                                         negate_zmm8r8(y));
	                                     }
	                                }	                                
	                                
	                        }    
	                        
	                        if(jint==2) {
	                           ei = negate_zmm8r8(ei);
	                        } 
	                        
	                        /*
	                           !
                                   !  To improve conditioning, rational approximations are expressed
                                   !  in terms of Chebyshev polynomials for 0 <= X < 6, and in
                                   !  continued fraction form for larger X.
                                   !
	                        */
	                        
	                }
	                else if(_mm512_cmp_pd_mask(x,six,_CMP_LT_OQ)) {
	                       
	                         t     = _mm512_add_pd(x,x);
	                         t     = _mm512_sub_pd(_mm512_div_pd(t,three),two);
	                         px[0] = zero;
	                         qx[0] = zero;
	                         px[1] = p[0];
	                         qx[1] = q[0];
	                         
	                         px[2] = _mm512_fmsub_pd(t,px[1],_mm512_add_pd(px[0],p[1]));
	                         qx[2] = _mm512_fmsub_pd(t,qx[1],_mm512_add_pd(qx[0],q[1]));
	                         px[3] = _mm512_fmsub_pd(t,px[2],_mm512_add_pd(px[1],p[2]));
	                         qx[3] = _mm512_fmsub_pd(t,qx[2],_mm512_add_pd(qx[1],q[2]));
	                         px[4] = _mm512_fmsub_pd(t,px[3],_mm512_add_pd(px[2],p[3]));
	                         qx[4] = _mm512_fmsub_pd(t,qx[3],_mm512_add_pd(qx[2],q[3]));
	                         px[5] = _mm512_fmsub_pd(t,px[4],_mm512_add_pd(px[3],p[4]));
	                         qx[5] = _mm512_fmsub_pd(t,qx[4],_mm512_add_pd(qx[3],q[4]));
	                         px[6] = _mm512_fmsub_pd(t,px[5],_mm512_add_pd(px[4],p[5]));
	                         qx[6] = _mm512_fmsub_pd(t,qx[5],_mm512_add_pd(qx[4],q[5]));
	                         px[7] = _mm512_fmsub_pd(t,px[6],_mm512_add_pd(px[5],p[6]));
	                         qx[7] = _mm512_fmsub_pd(t,qx[6],_mm512_add_pd(qx[5],q[6]));
	                         px[8] = _mm512_fmsub_pd(t,px[7],_mm512_add_pd(px[6],p[7]));
	                         qx[8] = _mm512_fmsub_pd(t,qx[7],_mm512_add_pd(qx[6],q[7]));
	                         __m512 tmp = _mm512_mul_pd(half,t);
	                         sump  = _mm512_fmsub_pd(tmp,px[9],_mm512_add_pd(px[8],p[9]));
	                         sumq  = _mm512_fmsub_pd(tmp,qx[9],_mm512_add_pd(qx[8],q[9]));
	                         frac  = _mm512_div_pd(sump,sumq);
	                         tmp    = _mm512_div_sub_pd(x,_mm512_div_pd(x01,x11));
	                         xmx0  = _mm512_sub_pd(tmp,x02);

                                 if(_mm512_cmp_pd_mask(p037,
                                                   _mm512_abs_pd(xmx0),_CMP_LE_OQ)) {
                                     
                                     __m512 tmp = _mm512_div_pd(x,x0);
                                     ei         = _mm512_fmadd_pd(frac,xmx0,xlog(tmp));
                                     if(jint==3) {
                                        ei = _mm512_mul_pd(xexp(negate_zmm8r8(x),ei));
                                     }                 
                                 }
                                 else {
                                    
                                      //Special approximation to ln(X/X0) for X close to X0. 
                                      
                                      y    = _mm512_div_pd(xmx0,_mm512_add_pd(x,x0));
                                      ysq  = _mm512_mul_pd(y,y);
                                      sump = plg[0]; 
                                      sumq = _mm512_add_pd(ysq,qlg[0]);
                                      
                                      sump = _mm512_fmadd_pd(sump,ysq,plg[1]);
                                      sumq = _mm512_fmadd_pd(sumq,ysq,qlg[1]);
                                      sump = _mm512_fmadd_pd(sump,ysq,plg[2]);
                                      sumq = _mm512_fmadd_pd(sumq,ysq,qlg[2]);
                                      sump = _mm512_fmadd_pd(sump,ysq,plg[3]);
                                      sumq = _mm512_fmadd_pd(sumq,ysq,qlg[3]);
                                      
                                      __m512 tmp = _mm512_fmadd_pd(sumq,
                                                               _mm512_add_pd(x,x0),frac);
                                      ei   = _mm512_mul_pd(_mm512_div_pd(sump,tmp),xmx0);
                                      
                                      if(jint==3) {
                                         ei = _mm512_mul_pd(xexp(negate_zmm8r8(x),ei));
                                      }
                                 }
	                }
	                else if(_mm512_cmp_pd_mask(x,twelve,_CMP_LT_OQ)) {
	                         
	                         frac = zero;
	                         
	                         frac = _mm512_div_pd(s[0],
	                                          _mm512_add_pd(r[0],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[1],
	                                          _mm512_add_pd(r[1],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[2],
	                                          _mm512_add_pd(r[2],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[3],
	                                          _mm512_add_pd(r[3],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[4],
	                                          _mm512_add_pd(r[4],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[5],
	                                          _mm512_add_pd(r[5],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[6],
	                                          _mm512_add_pd(r[6],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[7],
	                                          _mm512_add_pd(r[7],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(s[8],
	                                          _mm512_add_pd(r[8],
	                                                    _mm512_add_pd(frac)));
	                                                    
	                         ei   = _mm512_div_pd(_mm512_add_pd(r[9],frac),x);
	                         if(jint!=3) {
	                            ei = _mm512_mul_pd(ei,xexp(x));
	                         }
	                         
	                }
	                else if(_mm512_cmp_pd_mask(x,two4,_CMP_LE_OQ)) {
	                         
	                         frac = zero;
	                         
	                         frac = _mm512_div_pd(q1[0],
	                                          _mm512_add_pd(p1[0],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[1],
	                                          _mm512_add_pd(p1[1],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[2],
	                                          _mm512_add_pd(p1[2],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[3],
	                                          _mm512_add_pd(p1[3],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[4],
	                                          _mm512_add_pd(p1[4],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[5],
	                                          _mm512_add_pd(p1[5],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[6],
	                                          _mm512_add_pd(p1[6],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[7],
	                                          _mm512_add_pd(p1[7],
	                                                    _mm512_add_pd(frac)));
	                         frac = _mm512_div_pd(q1[8],
	                                          _mm512_add_pd(p1[8],
	                                                    _mm512_add_pd(frac)));
	                                                    
	                         ei   = _mm512_div_pd(_mm512_add_pd(p1[9],frac),x);
	                         if(jint!=3) {
	                            ei = _mm512_mul_pd(ei,xexp(x));
	                         }
	                         
	                }
	                else {
	                   
	                       if(_mm512_cmp_pd_mask(xmax,x,_CMP_LE_OQ) &&
	                          jint==3) {
	                          
	                           ei = xinf;  
	                       }
	                       else {
	                          
	                           y    = _mm512_div_pd(one,x);
	                           frac = zero;
	                           
	                           frac = _mm512_div_pd(q2[0],
	                                          _mm512_add_pd(p2[0],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[1],
	                                          _mm512_add_pd(p2[1],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[2],
	                                          _mm512_add_pd(p2[2],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[3],
	                                          _mm512_add_pd(p2[3],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[4],
	                                          _mm512_add_pd(p2[4],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[5],
	                                          _mm512_add_pd(p2[5],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[6],
	                                          _mm512_add_pd(p2[6],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[7],
	                                          _mm512_add_pd(p2[7],
	                                                    _mm512_add_pd(frac)));
	                           frac = _mm512_div_pd(q2[8],
	                                          _mm512_add_pd(p2[8],
	                                                    _mm512_add_pd(frac)));
	                                                    
	                           frac = _mm512_add_pd(p2[9],frac);
	                           ei   = _mm512_fmadd_pd(frac,
	                                              _mm512_mul_pd(y,y),y);
	                           if(jint!=3) {
	                              
	                              if(_mm512_cmp_pd_mask(x,
	                                       _mm512_sub_pd(xmax,two4),_CMP_LE_OQ)) {
	                                  
	                                  ei = _mm512_mul_pd(ei,xexp(x));    
	                                  // Calculation reformulated to avoid premature overflow.    
	                              }
	                              else {
	                                 
	                                  __m512 tmp = _mm512_sub_pd(x,fourty);
	                                  ei  = _mm512_mul_pd(_mm512_mul_pd(ei,xexp(tmp),exp40));
	                              }
	                           }
	                       }
	                }
	                 
	                 result = ei;
	                 return (result);
	                                
	         }
	         
	         
/*
    !*****************************************************************************80
!
!! CALCI0 computes various I0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the first kind
!    and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.
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
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
!    the argument must be less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = I0(x);
!    2, RESULT = exp(-x) * I0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, I0(x);
!    2, exp(-x) * I0(x);      
*/


                 
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calci0_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[15] = {_mm512_set1_pd(-5.2487866627945699800e-18),
	                                   _mm512_set1_pd(-1.5982226675653184646e-14),
                                           _mm512_set1_pd(-2.6843448573468483278e-11),
                                           _mm512_set1_pd(-3.0517226450451067446e-08), 
                                           _mm512_set1_pd(-2.5172644670688975051e-05),
                                           _mm512_set1_pd(-1.5453977791786851041e-02), 
                                           _mm512_set1_pd(-7.0935347449210549190e+00),
                                           _mm512_set1_pd(-2.4125195876041896775e+03), 
                                           _mm512_set1_pd(-5.9545626019847898221e+05),
                                           _mm512_set1_pd(-1.0313066708737980747e+08), 
                                           _mm512_set1_pd(-1.1912746104985237192e+10),
                                           _mm512_set1_pd(-8.4925101247114157499e+11), 
                                           _mm512_set1_pd(-3.2940087627407749166e+13),
                                           _mm512_set1_pd(-5.5050369673018427753e+14), 
                                           _mm512_set1_pd(-2.2335582639474375249e+15)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[5] =  {_mm512_set1_pd(3.7277560179962773046e+03), 
	                                   _mm512_set1_pd(6.5158506418655165707e+06), 
                                           _mm512_set1_pd(-6.5626560740833869295e+09), 
                                           _mm512_set1_pd(3.7604188704092954661e+12), 
                                           _mm512_set1_pd(-9.7087946179594019126e+14)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[8] = {_mm512_set1_pd(3.9843750000000000000e-01), 
	                                   _mm512_set1_pd(2.9205384596336793945e+00), 
                                           _mm512_set1_pd(-2.4708469169133954315e+00), 
                                           _mm512_set1_pd(4.7914889422856814203e-01), 
                                           _mm512_set1_pd(-3.7384991926068969150e-03),
                                           _mm512_set1_pd(-2.6801520353328635310e-03), 
                                           _mm512_set1_pd(9.9168777670983678974e-05),
                                           _mm512_set1_pd(-2.1877128189032726730e-06)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d qq[7] = {_mm512_set1_pd(3.1446690275135491500e+01), 
	                                   _mm512_set1_pd(8.5539563258012929600e+01), 
                                           _mm512_set1_pd(-6.0228002066743340583e+01), 
                                           _mm512_set1_pd(1.3982595353892851542e+01), 
                                           _mm512_set1_pd(-1.1151759188741312645e+00), 
                                           _mm512_set1_pd(3.2547697594819615062e-02), 
                                           _mm512_set1_pd(-5.5194330231005480228e-04)};
                         const __m512d one   = _mm512_set1_pd(1.0e+0);
                         const __m512d one5  = _mm512_set1_pd(15.0e+0);  
                         const __m512d exp40 = _mm512_set1_pd(2.353852668370199854e+17);
                         const __m512d forty = _mm512_set1_pd(40.0e+0);
                         const __m512d rec15 = _mm512_set1_pd(6.6666666666666666666e-2);
                         const __m512d two25 = _mm512_set1_pd(225.0e+0);
                         const __m512d xsmall= _mm512_set1_pd(5.55e-17);
                         const __m512d xinf  = _mm512_set1_pd(1.79e+308);
                         const __m512d xmax  = _mm512_set1_pd(713.986e+0);
                         register __m512d a,b,sump,sumq,x,xx;
                         register __m512d result;
                         
                         x = _mm512_abs_pd(arg);
                         if(_mm512_cmp_pd_mask(x,xsmall,_CMP_LT_OQ)) {
                            
                             result = one;
                             /*
                                  XSMALL <= ABS(ARG) < 15.0.
                             */
                         }   
                         else if(_mm512_cmp_pd_mask(x,one5,_CMP_LT_OQ)) {
                              
                              xx   = _mm512_mul_pd(x,x);
                              sump = p[0];
                              sump = _mm512_fmadd_pd(sump,xx,p[1]);
                              sump = _mm512_fmadd_pd(sump,xx,p[2]);
                              sump = _mm512_fmadd_pd(sump,xx,p[3]);
                              sump = _mm512_fmadd_pd(sump,xx,p[4]);
                              sump = _mm512_fmadd_pd(sump,xx,p[5]);
                              sump = _mm512_fmadd_pd(sump,xx,p[6]);
                              sump = _mm512_fmadd_pd(sump,xx,p[7]);
                              sump = _mm512_fmadd_pd(sump,xx,p[8]);
                              sump = _mm512_fmadd_pd(sump,xx,p[9]);
                              sump = _mm512_fmadd_pd(sump,xx,p[10]);
                              sump = _mm512_fmadd_pd(sump,xx,p[11]);
                              sump = _mm512_fmadd_pd(sump,xx,p[12]);
                              sump = _mm512_fmadd_pd(sump,xx,p[13]);
                              sump = _mm512_fmadd_pd(sump,xx,p[14]);
                              xx   = _mm512_sub_pd(xx,two25);
                              
                              /*register __m512d t0 = _mm512_add_pd(xx,q[0]);
                              register __m512d t1 = _mm512_add_pd(xx,q[1]);
                              register __m512d t2 = _mm512_add_pd(xx,q[2]);
                              register __m512d t3 = _mm512_add_pd(xx,q[3]);
                              register __m512d t4 = _mm512_add_pd(xx,q[4]);
                              sumq                = _mm512_mul_pd(
                                                           _mm512_mul_pd(
                                                                   _mm512_mul_pd(t0,t1),
                                                                          _mm512_mul_pd(t2,t3)),t4);*/
                              sumq    = _mm512_fmadd_pd(
                                                    _mm512_add_pd(xx,q[0]),xx,
                                                               _mm512_fmadd_pd(q[1],xx,
                                                                      _mm512_fmadd_pd(q[2],xx,
                                                                             _mm512_fmadd_pd(q[3],xx,
                                                                                 _mm512_fmadd_pd(q[4],xx,q[5])))));
                                                                                                                                
                              result              = _mm512_div_pd(sump,sumq);
                              
                              if(jint==2) {
                                  result = _mm512_mul_pd(result,
                                                    xexp(negate_zmm8r8(x)));
                              }
                              
                         }    
                         else if(_mm512_cmp_pd_mask(one5,x,_CMP_LE_OQ)) {
                                 
                                 if(jint==1 && 
                                     _mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                      result = xinf;    
                                 }
                                 else {
                                    
                                     xx = _mm512_div_pd(one,
                                                    _mm512_sub_pd(x,rec15));
                                     sump = _mm512_fmadd_pd(pp[0],xx,
                                                   _mm512_fmadd_pd(pp[1],xx,
                                                         _mm512_fmadd_pd(pp[2],xx,
                                                               _mm512_fmadd_pd(pp[3],xx,
                                                                     _mm512_fmadd_pd(pp[4],xx,
                                                                           _mm512_fmadd_pd(pp[5],xx,
                                                                                 _mm512_fmadd_pd(pp[6],xx,
                                                                                       _mm512_mul_pd(pp[7],xx))))))));
                                     sumq = _mm512_fmadd_pd(
                                                  _mm512_add_pd(xx,q[0]),xx,
                                                       _mm512_fmadd_pd(q[1],xx,
                                                             _mm512_fmadd_pd(q[2],xx,
                                                                  _mm512_fmadd_pd(q[3],xx,
                                                                        _mm512_fmadd_pd(q[4],xx,
                                                                            _mm512_fmadd_pd(q[5],xx,q[6]))))));
                                     result = _mm512_div_pd(sump,sumq);
                                     
                                     if(jint==2) {
                                        register __m512d tmp = _mm512_sqrt_pd(x);
                                        result  = _mm512_div_pd(_mm512_sub_pd(result,pp[0]),tmp);
                                     } 
                                     else { 
                                         
                                          if(_mm512_cmp_pd_mask(x,
                                                    _mm512_sub_pd(xmax,one5),_CMP_LE_OQ)) {
                                             a = xexp(x);
                                             b = one;             
                                          }
                                          else {
                                             a = xexp(_mm512_sub_pd(x,forty));
                                             b = exp40;
                                          }
                                          
                                          register __m512 tmp = _mm512_sqrt_pd(x);
                                          result  = _mm512_mul_pd(_mm512_fmsub_pd(result,a,
                                                                          _mm512_mul_pd(pp[0],a)),a);
                                          result  = _mm512_mul_pd(_mm512_div_pd(result,tmp),b);
                                     }
                                     
                                     /*
                                         Calculation reformulated to avoid premature overflow.
                                     */                                            
                                     
                                 }
                         }
                         
                         return (result);
	         }   
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calci0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[15] = {_mm512_set1_pd(-5.2487866627945699800e-18),
	                                   _mm512_set1_pd(-1.5982226675653184646e-14),
                                           _mm512_set1_pd(-2.6843448573468483278e-11),
                                           _mm512_set1_pd(-3.0517226450451067446e-08), 
                                           _mm512_set1_pd(-2.5172644670688975051e-05),
                                           _mm512_set1_pd(-1.5453977791786851041e-02), 
                                           _mm512_set1_pd(-7.0935347449210549190e+00),
                                           _mm512_set1_pd(-2.4125195876041896775e+03), 
                                           _mm512_set1_pd(-5.9545626019847898221e+05),
                                           _mm512_set1_pd(-1.0313066708737980747e+08), 
                                           _mm512_set1_pd(-1.1912746104985237192e+10),
                                           _mm512_set1_pd(-8.4925101247114157499e+11), 
                                           _mm512_set1_pd(-3.2940087627407749166e+13),
                                           _mm512_set1_pd(-5.5050369673018427753e+14), 
                                           _mm512_set1_pd(-2.2335582639474375249e+15)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[5] =  {_mm512_set1_pd(3.7277560179962773046e+03), 
	                                   _mm512_set1_pd(6.5158506418655165707e+06), 
                                           _mm512_set1_pd(-6.5626560740833869295e+09), 
                                           _mm512_set1_pd(3.7604188704092954661e+12), 
                                           _mm512_set1_pd(-9.7087946179594019126e+14)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[8] = {_mm512_set1_pd(3.9843750000000000000e-01), 
	                                   _mm512_set1_pd(2.9205384596336793945e+00), 
                                           _mm512_set1_pd(-2.4708469169133954315e+00), 
                                           _mm512_set1_pd(4.7914889422856814203e-01), 
                                           _mm512_set1_pd(-3.7384991926068969150e-03),
                                           _mm512_set1_pd(-2.6801520353328635310e-03), 
                                           _mm512_set1_pd(9.9168777670983678974e-05),
                                           _mm512_set1_pd(-2.1877128189032726730e-06)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d qq[7] = {_mm512_set1_pd(3.1446690275135491500e+01), 
	                                   _mm512_set1_pd(8.5539563258012929600e+01), 
                                           _mm512_set1_pd(-6.0228002066743340583e+01), 
                                           _mm512_set1_pd(1.3982595353892851542e+01), 
                                           _mm512_set1_pd(-1.1151759188741312645e+00), 
                                           _mm512_set1_pd(3.2547697594819615062e-02), 
                                           _mm512_set1_pd(-5.5194330231005480228e-04)};
                         const __m512d one   = _mm512_set1_pd(1.0e+0);
                         const __m512d one5  = _mm512_set1_pd(15.0e+0);  
                         const __m512d exp40 = _mm512_set1_pd(2.353852668370199854e+17);
                         const __m512d forty = _mm512_set1_pd(40.0e+0);
                         const __m512d rec15 = _mm512_set1_pd(6.6666666666666666666e-2);
                         const __m512d two25 = _mm512_set1_pd(225.0e+0);
                         const __m512d xsmall= _mm512_set1_pd(5.55e-17);
                         const __m512d xinf  = _mm512_set1_pd(1.79e+308);
                         const __m512d xmax  = _mm512_set1_pd(713.986e+0);
                         register __m512d arg,a,b,sump,sumq,x,xx;
                         register __m512d result;
                         
                         arg = _mm512_load_pd(&parg[0]);
                         x = _mm512_abs_pd(arg);
                         if(_mm512_cmp_pd_mask(x,xsmall,_CMP_LT_OQ)) {
                            
                             result = one;
                             /*
                                  XSMALL <= ABS(ARG) < 15.0.
                             */
                         }   
                         else if(_mm512_cmp_pd_mask(x,one5,_CMP_LT_OQ)) {
                              
                              xx   = _mm512_mul_pd(x,x);
                              sump = p[0];
                              sump = _mm512_fmadd_pd(sump,xx,p[1]);
                              sump = _mm512_fmadd_pd(sump,xx,p[2]);
                              sump = _mm512_fmadd_pd(sump,xx,p[3]);
                              sump = _mm512_fmadd_pd(sump,xx,p[4]);
                              sump = _mm512_fmadd_pd(sump,xx,p[5]);
                              sump = _mm512_fmadd_pd(sump,xx,p[6]);
                              sump = _mm512_fmadd_pd(sump,xx,p[7]);
                              sump = _mm512_fmadd_pd(sump,xx,p[8]);
                              sump = _mm512_fmadd_pd(sump,xx,p[9]);
                              sump = _mm512_fmadd_pd(sump,xx,p[10]);
                              sump = _mm512_fmadd_pd(sump,xx,p[11]);
                              sump = _mm512_fmadd_pd(sump,xx,p[12]);
                              sump = _mm512_fmadd_pd(sump,xx,p[13]);
                              sump = _mm512_fmadd_pd(sump,xx,p[14]);
                              xx   = _mm512_sub_pd(xx,two25);
                              
                              /*register __m512d t0 = _mm512_add_pd(xx,q[0]);
                              register __m512d t1 = _mm512_add_pd(xx,q[1]);
                              register __m512d t2 = _mm512_add_pd(xx,q[2]);
                              register __m512d t3 = _mm512_add_pd(xx,q[3]);
                              register __m512d t4 = _mm512_add_pd(xx,q[4]);
                              sumq                = _mm512_mul_pd(
                                                           _mm512_mul_pd(
                                                                   _mm512_mul_pd(t0,t1),
                                                                          _mm512_mul_pd(t2,t3)),t4);*/
                              sumq    = _mm512_fmadd_pd(
                                                    _mm512_add_pd(xx,q[0]),xx,
                                                               _mm512_fmadd_pd(q[1],xx,
                                                                      _mm512_fmadd_pd(q[2],xx,
                                                                             _mm512_fmadd_pd(q[3],xx,q[4]))));
                                                                                 
                                                                                                                                
                              result              = _mm512_div_pd(sump,sumq);
                              
                              if(jint==2) {
                                  result = _mm512_mul_pd(result,
                                                    xexp(negate_zmm8r8(x)));
                              }
                              
                         }    
                         else if(_mm512_cmp_pd_mask(one5,x,_CMP_LE_OQ)) {
                                 
                                 if(jint==1 && 
                                     _mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                      result = xinf;    
                                 }
                                 else {
                                    
                                     xx = _mm512_div_pd(one,
                                                    _mm512_sub_pd(x,rec15));
                                     sump = _mm512_fmadd_pd(pp[0],xx,
                                                   _mm512_fmadd_pd(pp[1],xx,
                                                         _mm512_fmadd_pd(pp[2],xx,
                                                               _mm512_fmadd_pd(pp[3],xx,
                                                                     _mm512_fmadd_pd(pp[4],xx,
                                                                           _mm512_fmadd_pd(pp[5],xx,
                                                                                 _mm512_fmadd_pd(pp[6],xx,
                                                                                       _mm512_mul_pd(pp[7],xx))))))));
                                     sumq = _mm512_fmadd_pd(
                                                  _mm512_add_pd(xx,qq[0]),xx,
                                                       _mm512_fmadd_pd(qq[1],xx,
                                                             _mm512_fmadd_pd(qq[2],xx,
                                                                  _mm512_fmadd_pd(qq[3],xx,
                                                                        _mm512_fmadd_pd(qq[4],xx,
                                                                            _mm512_fmadd_pd(qq[5],xx,qq[6]))))));
                                     result = _mm512_div_pd(sump,sumq);
                                     
                                     if(jint==2) {
                                        register __m512d tmp = _mm512_sqrt_pd(x);
                                        result  = _mm512_div_pd(_mm512_sub_pd(result,pp[0]),tmp);
                                     } 
                                     else { 
                                         
                                          if(_mm512_cmp_pd_mask(x,
                                                    _mm512_sub_pd(xmax,one5),_CMP_LE_OQ)) {
                                             a = xexp(x);
                                             b = one;             
                                          }
                                          else {
                                             a = xexp(_mm512_sub_pd(x,forty));
                                             b = exp40;
                                          }
                                          
                                          register __m512 tmp = _mm512_sqrt_pd(x);
                                          result  = _mm512_mul_pd(_mm512_fmsub_pd(result,a,
                                                                          _mm512_mul_pd(pp[0],a)),a);
                                          result  = _mm512_mul_pd(_mm512_div_pd(result,tmp),b);
                                     }
                                     
                                     /*
                                         Calculation reformulated to avoid premature overflow.
                                     */                                            
                                     
                                 }
                         }
                         
                         return (result);
	         }   
	         
	         
/*
!*****************************************************************************80
!
!! CALCI1 computes various I1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functioons of the first kind
!    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.
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
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
!    the argument must be less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = I1(x);
!    2, RESULT = exp(-x) * I1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, I1(x);
!    2, exp(-x) * I1(x);
*/	


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calci1_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {    
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[15] = {_mm512_set1_pd(-1.9705291802535139930e-19),
	                                   _mm512_set1_pd(-6.5245515583151902910e-16), 
                                           _mm512_set1_pd(-1.1928788903603238754e-12),
                                           _mm512_set1_pd(-1.4831904935994647675e-09), 
                                           _mm512_set1_pd(-1.3466829827635152875e-06),
                                           _mm512_set1_pd(-9.1746443287817501309e-04), 
                                           _mm512_set1_pd(-4.7207090827310162436e-01),
                                           _mm512_set1_pd(-1.8225946631657315931e+02), 
                                           _mm512_set1_pd(-5.1894091982308017540e+04),
                                           _mm512_set1_pd(-1.0588550724769347106e+07), 
                                           _mm512_set1_pd(-1.4828267606612366099e+09),
                                           _mm512_set1_pd(-1.3357437682275493024e+11), 
                                           _mm512_set1_pd(-6.9876779648010090070e+12),
                                           _mm512_set1_pd(-1.7732037840791591320e+14), 
                                           _mm512_set1_pd(-1.4577180278143463643e+15)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[5] =   {_mm512_set1_pd(-4.0076864679904189921e+03), 
	                                    _mm512_set1_pd(7.4810580356655069138e+06), 
                                            _mm512_set1_pd(-8.0059518998619764991e+09), 
                                            _mm512_set1_pd(4.8544714258273622913e+12), 
                                            _mm512_set1_pd(-1.3218168307321442305e+15)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[8] =  {_mm512_set1_pd(6.0437159056137600000e-02), 
	                                    _mm512_set1_pd(4.5748122901933459000e-01),
                                            _mm512_set1_pd(-4.2843766903304806403e-01), 
                                            _mm512_set1_pd(9.7356000150886612134e-02), 
                                            _mm512_set1_pd(-3.2457723974465568321e-03),
                                            _mm512_set1_pd(-3.6395264712121795296e-04), 
                                            _mm512_set1_pd(1.6258661867440836395e-05),
                                            _mm512_set1_pd(-3.6347578404608223492e-07)};
                           __ATTR_ALIGN__(64) const static 
	                  __m512d qq[6] =  {_mm512_set1_pd(-3.8806586721556593450e+00), 
	                                    _mm512_set1_pd(3.2593714889036996297e+00), 
                                            _mm512_set1_pd(-8.5017476463217924408e-01), 
                                            _mm512_set1_pd(7.4212010813186530069e-02), 
                                            _mm512_set1_pd(-2.2835624489492512649e-03), 
                                            _mm512_set1_pd(3.7510433111922824643e-05)};
                          const __m512d one   = _mm512_set1_pd(1.0e+0);
                          const __m512d one5  = _mm512_set1_pd(15.0e+0);
                          const __m512d exp40 = _mm512_set1_pd(2.353852668370199854e+17);
                          const __m512d forty = _mm512_set1_pd(40.0e+0);
                          const __m512d rec15 = _mm512_set1_pd(6.6666666666666666666e-2);
                          const __m512d two25 = _mm512_set1_pd(225.0e+0);
                          const __m512d half  = _mm512_set1_pd(0.5e+0);
                          const __m512d zero  = _mm512_set1_pd(0.0e+0);
//!
//!  Machine-dependent constants
//!
                          const __m512d xsmall = _mm512_set1_pd(5.55e-17);
                          const __m512d xinf   = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax   = _mm512_set1_pd(713.987e+0);
                          const __m512d pbar   = _mm512_set1_pd(3.98437500e-01);  
                          register __m512d a,b,result,sump,sumq;
                          register __m512d x,xx;
                          
                          x = _mm512_abs_pd(arg);
                          //
                          // Return for ABS(ARG) < XSMALL.
                          //   
                          if(_mm512_cmp_pd_mask(x,xsmall,_CMP_LT_OQ)) {
                          
                               result = _mm512_mul_pd(half,x);
                          }    
                          else if(_mm512_cmp_pd_mask(x,one5,_CMP_LT_OQ)) {
                          
                               xx    = _mm512_mul_pd(x,x);
                               sump  = p[0];
                               sump  = _mm512_fmadd_pd(sump,xx,p[1]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[2]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[3]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[4]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[5]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[6]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[7]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[8]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[9]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[10]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[11]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[12]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[13]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[14]);
                               xx    = _mm512_sub_pd(xx,two25);
                               
                               sumq  = _mm512_fmadd_pd(
                                             _mm512_add_pd(xx,q[0]),xx,
                                                   _mm512_fmadd_pd(q[1],xx,
                                                         _mm512_fmadd_pd(q[2],xx,
                                                               _mm512_fmadd_pd(q[3],xx,q[4]))));
                               result = _mm512_mul_pd(_mm512_div_pd(sump,sumq),x);
                               
                               if(jint==2) {
                                  result = _mm512_mul_pd(result,xexp(negate_zmm8r8(x));
                               }
                     
                         }   
                         else if(jint==1 &&
                                 _mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                 
                                 result = xinf;        
                        }                         
                        else {
                                
                               xx = _mm512_sub_pd(_mm512_div_pd(one,x),rec15);
                               sump = _mm512_fmadd_pd(pp[0],xx,
                                            _mm512_fmadd_pd(pp[1],xx,
                                                  _mm512_fmadd_pd(pp[2],xx,
                                                        _mm512_fmadd_pd(pp[3],xx,
                                                              _mm512_fmadd_pd(pp[4],xx,
                                                                    _mm512_fmadd_pd(pp[5],xx,
                                                                          _mm512_fmadd_pd(pp[6],xx,
                                                                                _mm512_mul_pd(pp[7],xx))))))));
                               sumq = _mm512_fmadd_pd(
                                             _mm512_add_pd(xx,qq[0]),xx,
                                                   _mm512_fmadd_pd(qq[1],xx,
                                                         _mm512_fmadd_pd(qq[2],xx,
                                                               _mm512_fmadd_pd(qq[3],xx,
                                                                     _mm512_fmadd_pd(qq[4],xx,qq[5])))));
                               result = _mm512_div_pd(sump,sumq);
                               if(jint!=1) {
                                  register __m512d tmp = _mm512_sqrt_pd(x);
                                  result = _mm512_div_pd(_mm512_add_pd(result,pbar),tmp);
                               }  
                               else {
                                  const __mmask16 m = 
                                             _mm512_cmp_pd_mask(
                                                         _mm512_sub_pd(xmax,one5),_CMP_LT_OQ);
                                  a  = _mm512_mask_blend_pd(m,xexp(x),
                                                       xexp(_mm512_sub_pd(x,forty)))
                                  b  = _mm512_mask_blend_pd(m,one,exp40);
                                  register __m512d tmp = _mm512_fmadd_pd(result,a,
                                                                     _mm512_mul_pd(pbar,a));
                                  result               = _mm512_mul_pd(_mm512_div_pd(tmp,
                                                                          _mm512_sqrt_pd(x)),b);                        
                                                            
                               }                                         
                        } 
                        
                        if(_mm512_cmp_pd_mask(arg,zero,_CMP_LT_OQ)) {
                            
                            result = negate_zmm8r8(result);
                        }  
                        
                        return (result);     
	          }   
	          
	          
	              
	          
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calci1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                 const int32_t jint) {    
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[15] = {_mm512_set1_pd(-1.9705291802535139930e-19),
	                                   _mm512_set1_pd(-6.5245515583151902910e-16), 
                                           _mm512_set1_pd(-1.1928788903603238754e-12),
                                           _mm512_set1_pd(-1.4831904935994647675e-09), 
                                           _mm512_set1_pd(-1.3466829827635152875e-06),
                                           _mm512_set1_pd(-9.1746443287817501309e-04), 
                                           _mm512_set1_pd(-4.7207090827310162436e-01),
                                           _mm512_set1_pd(-1.8225946631657315931e+02), 
                                           _mm512_set1_pd(-5.1894091982308017540e+04),
                                           _mm512_set1_pd(-1.0588550724769347106e+07), 
                                           _mm512_set1_pd(-1.4828267606612366099e+09),
                                           _mm512_set1_pd(-1.3357437682275493024e+11), 
                                           _mm512_set1_pd(-6.9876779648010090070e+12),
                                           _mm512_set1_pd(-1.7732037840791591320e+14), 
                                           _mm512_set1_pd(-1.4577180278143463643e+15)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[5] =   {_mm512_set1_pd(-4.0076864679904189921e+03), 
	                                    _mm512_set1_pd(7.4810580356655069138e+06), 
                                            _mm512_set1_pd(-8.0059518998619764991e+09), 
                                            _mm512_set1_pd(4.8544714258273622913e+12), 
                                            _mm512_set1_pd(-1.3218168307321442305e+15)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[8] =  {_mm512_set1_pd(6.0437159056137600000e-02), 
	                                    _mm512_set1_pd(4.5748122901933459000e-01),
                                            _mm512_set1_pd(-4.2843766903304806403e-01), 
                                            _mm512_set1_pd(9.7356000150886612134e-02), 
                                            _mm512_set1_pd(-3.2457723974465568321e-03),
                                            _mm512_set1_pd(-3.6395264712121795296e-04), 
                                            _mm512_set1_pd(1.6258661867440836395e-05),
                                            _mm512_set1_pd(-3.6347578404608223492e-07)};
                           __ATTR_ALIGN__(64) const static 
	                  __m512d qq[6] =  {_mm512_set1_pd(-3.8806586721556593450e+00), 
	                                    _mm512_set1_pd(3.2593714889036996297e+00), 
                                            _mm512_set1_pd(-8.5017476463217924408e-01), 
                                            _mm512_set1_pd(7.4212010813186530069e-02), 
                                            _mm512_set1_pd(-2.2835624489492512649e-03), 
                                            _mm512_set1_pd(3.7510433111922824643e-05)};
                          const __m512d one   = _mm512_set1_pd(1.0e+0);
                          const __m512d one5  = _mm512_set1_pd(15.0e+0);
                          const __m512d exp40 = _mm512_set1_pd(2.353852668370199854e+17);
                          const __m512d forty = _mm512_set1_pd(40.0e+0);
                          const __m512d rec15 = _mm512_set1_pd(6.6666666666666666666e-2);
                          const __m512d two25 = _mm512_set1_pd(225.0e+0);
                          const __m512d half  = _mm512_set1_pd(0.5e+0);
                          const __m512d zero  = _mm512_set1_pd(0.0e+0);
//!
//!  Machine-dependent constants
//!
                          const __m512d xsmall = _mm512_set1_pd(5.55e-17);
                          const __m512d xinf   = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax   = _mm512_set1_pd(713.987e+0);
                          const __m512d pbar   = _mm512_set1_pd(3.98437500e-01);  
                          register __m512d arg = _mm512_load_pd(&parg[0]);
                          register __m512d a,b,result,sump,sumq;
                          register __m512d x,xx;
                          
                          x = _mm512_abs_pd(arg);
                          //
                          // Return for ABS(ARG) < XSMALL.
                          //   
                          if(_mm512_cmp_pd_mask(x,xsmall,_CMP_LT_OQ)) {
                          
                               result = _mm512_mul_pd(half,x);
                          }    
                          else if(_mm512_cmp_pd_mask(x,one5,_CMP_LT_OQ)) {
                          
                               xx    = _mm512_mul_pd(x,x);
                               sump  = p[0];
                               sump  = _mm512_fmadd_pd(sump,xx,p[1]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[2]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[3]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[4]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[5]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[6]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[7]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[8]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[9]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[10]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[11]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[12]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[13]);
                               sump  = _mm512_fmadd_pd(sump,xx,p[14]);
                               xx    = _mm512_sub_pd(xx,two25);
                               
                               sumq  = _mm512_fmadd_pd(
                                             _mm512_add_pd(xx,q[0]),xx,
                                                   _mm512_fmadd_pd(q[1],xx,
                                                         _mm512_fmadd_pd(q[2],xx,
                                                               _mm512_fmadd_pd(q[3],xx,q[4]))));
                               result = _mm512_mul_pd(_mm512_div_pd(sump,sumq),x);
                               
                               if(jint==2) {
                                  result = _mm512_mul_pd(result,xexp(negate_zmm8r8(x));
                               }
                     
                         }   
                         else if(jint==1 &&
                                 _mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                 
                                 result = xinf;        
                        }                         
                        else {
                                
                               xx = _mm512_sub_pd(_mm512_div_pd(one,x),rec15);
                               sump = _mm512_fmadd_pd(pp[0],xx,
                                            _mm512_fmadd_pd(pp[1],xx,
                                                  _mm512_fmadd_pd(pp[2],xx,
                                                        _mm512_fmadd_pd(pp[3],xx,
                                                              _mm512_fmadd_pd(pp[4],xx,
                                                                    _mm512_fmadd_pd(pp[5],xx,
                                                                          _mm512_fmadd_pd(pp[6],xx,
                                                                                _mm512_mul_pd(pp[7],xx))))))));
                               sumq = _mm512_fmadd_pd(
                                             _mm512_add_pd(xx,qq[0]),xx,
                                                   _mm512_fmadd_pd(qq[1],xx,
                                                         _mm512_fmadd_pd(qq[2],xx,
                                                               _mm512_fmadd_pd(qq[3],xx,
                                                                     _mm512_fmadd_pd(qq[4],xx,qq[5])))));
                               result = _mm512_div_pd(sump,sumq);
                               if(jint!=1) {
                                  register __m512d tmp = _mm512_sqrt_pd(x);
                                  result = _mm512_div_pd(_mm512_add_pd(result,pbar),tmp);
                               }  
                               else {
                                  const __mmask16 m = 
                                             _mm512_cmp_pd_mask(
                                                         _mm512_sub_pd(xmax,one5),_CMP_LT_OQ);
                                  a  = _mm512_mask_blend_pd(m,xexp(x),
                                                       xexp(_mm512_sub_pd(x,forty)))
                                  b  = _mm512_mask_blend_pd(m,one,exp40);
                                  register __m512d tmp = _mm512_fmadd_pd(result,a,
                                                                     _mm512_mul_pd(pbar,a));
                                  result               = _mm512_mul_pd(_mm512_div_pd(tmp,
                                                                          _mm512_sqrt_pd(x)),b);                        
                                                            
                               }                                         
                        } 
                        
                        if(_mm512_cmp_pd_mask(arg,zero,_CMP_LT_OQ)) {
                            
                            result = negate_zmm8r8(result);
                        }  
                        
                        return (result);     
	          }   
	          
	        
/*
*****************************************************************************80
!
!! CALCK0 computes various K0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
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
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K0(x);
!    2, RESULT = exp(x) * K0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
*/
	               
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calck0_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[6]  = {_mm512_set1_pd(5.8599221412826100000e-04), 
	                                  _mm512_set1_pd(1.3166052564989571850e-01), 
                                          _mm512_set1_pd(1.1999463724910714109e+01), 
                                          _mm512_set1_pd(4.6850901201934832188e+02), 
                                          _mm512_set1_pd(5.9169059852270512312e+03), 
                                          _mm512_set1_pd(2.4708152720399552679e+03)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[2]   = {_mm512_set1_pd(-2.4994418972832303646e+02), 
	                                  _mm512_set1_pd(2.1312714303849120380e+04)};
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d f[4]   = {_mm512_set1_pd(-1.6414452837299064100e+00),
	                                  _mm512_set1_pd(-2.9601657892958843866e+02), 
                                          _mm512_set1_pd(-1.7733784684952985886e+04),
                                          _mm512_set1_pd(-4.0320340761145482298e+05)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d g[3]   = {_mm512_set1_pd(-2.5064972445877992730e+02), 
	                                  _mm512_set1_pd(2.9865713163054025489e+04), 
                                          _mm512_set1_pd(-1.6128136304458193998e+06)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[10] = {_mm512_set1_pd(1.1394980557384778174e+02), 
	                                    _mm512_set1_pd(3.6832589957340267940e+03), 
                                            _mm512_set1_pd(3.1075408980684392399e+04), 
                                            _mm512_set1_pd(1.0577068948034021957e+05), 
                                            _mm512_set1_pd(1.7398867902565686251e+05), 
                                            _mm512_set1_pd(1.5097646353289914539e+05), 
                                            _mm512_set1_pd(7.1557062783764037541e+04), 
                                            _mm512_set1_pd(1.8321525870183537725e+04), 
                                            _mm512_set1_pd(2.3444738764199315021e+03), 
                                            _mm512_set1_pd(1.1600249425076035558e+02)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d qq[10] = {_mm512_set1_pd(2.0013443064949242491e+02), 
	                                    _mm512_set1_pd(4.4329628889746408858e+03), 
                                            _mm512_set1_pd(3.1474655750295278825e+04), 
                                            _mm512_set1_pd(9.7418829762268075784e+04), 
                                            _mm512_set1_pd(1.5144644673520157801e+05), 
                                            _mm512_set1_pd(1.2689839587977598727e+05), 
                                            _mm512_set1_pd(5.8824616785857027752e+04), 
                                            _mm512_set1_pd(1.4847228371802360957e+04), 
                                            _mm512_set1_pd(1.8821890840982713696e+03), 
                                            _mm512_set1_pd(9.2556599177304839811e+01)};

                          const __m512d one    = _mm512_set1_pd(1.0e+0);
                          const __m512d zero   = _mm512_set1_pd(0.0e+0);
                          const __m512d xsmall = _mm512_set1_pd(1.11e-16);
                          const __m512d xinf   = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax   = _mm512_set1_pd(705.342e+0);
                          register __m512d result,sumf,sumg,sump,sumq,temp;
                          register __m512d x,xx;
                          
                          x = arg;
                          if(_mm512_cmp_pd_mask(zero,x,_CMP_LT_OQ)) {
                             
                               if(_mm512_cmp_pd_mask(x,one,_CMP_LE_OQ)) {
                                   
                                   temp = xlog(x);
                                   
                                   if(_mm512_cmp_pd_mask(x,small,_CMP_LT_OQ)) {
                                       
                                       result = _mm512_sub_pd(_mm512_div_pd(p[5],p[1]),temp);
                                   }
                                   else {
                                       
                                        xx   = _mm512_mul_pd(x,x);
                                        sump = _mm512_fmadd_pd(pp[0],xx,
                                                      _mm512_fmadd_pd(pp[1],xx,
                                                            _mm512_fmadd_pd(pp[2],xx,
                                                                  _mm512_fmadd_pd(pp[3],xx,
                                                                        _mm512_fmadd_pd(pp[4],xx,p[5])))));
                                        sumq = _mm512_fmadd_pd(xx,q[0],
                                                      _mm512_add_pd(xx,q[1]));
                                        sumf = _mm512_fmadd_pd(f[0],xx,
                                                      _mm512_fmadd_pd(f[1],xx,
                                                             _mm512_fmadd_pd(f[2],xx,f[3])));
                                        sumg = _mm512_fmadd_pd(
                                                       _mm512_add_pd(xx,g[0]),xx,
                                                                 _mm512_fmadd_pd(g[1],xx,g[2]));
                                        register __m512d t0 = _mm512_div_pd(sump,sumq);
                                        register __m512d t1 = _mm512_sub_pd(_mm512_div_pd(temp,sumg),temp);
                                        result   = _mm512_sub_pd(t0,
                                                           _mm512_mul_pd(xx,
                                                                  _mm512_mul_pd(sumf,t1)));
                                        if(jint==2) {
                                            result = _mm512_mul_pd(result,xexp(x));
                                        }
                                   }
                               }
                               else if(jint==1 && 
                                       _mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                     
                                       result = zero;       
                               }
                               else {
                                       
                                      xx   = _mm512_div_pd(one,x);
                                      sumq = xx;
                                      sump = p[0];
                                      sump = _mm512_fmadd_pd(sump,xx,p[1]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[0]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[2]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[1]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[3]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[2]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[4]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[3]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[5]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[4]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[6]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[5]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[7]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[6]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[8]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[7]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[9]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[8]));
                                      sumq = _mm512_add_pd(sumq,qq[9]);
                                      register __m512d t0 = _mm512_sqrt_pd(x);
                                      register __m512d t1 = _mm512_div_pd(sump,sumq);
                                      result = _mm512_div_pd(t1,t0);
                                      if(jint==1) {
                                         result = _mm512_mul_pd(result,xexp(negate_zmm8r8(x)));
                                      }
                               }
                          }  
                          else {
                                  result = xinf;
                          }  
                          
                          return (result);                  
	        }    
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calck0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[6]  = {_mm512_set1_pd(5.8599221412826100000e-04), 
	                                  _mm512_set1_pd(1.3166052564989571850e-01), 
                                          _mm512_set1_pd(1.1999463724910714109e+01), 
                                          _mm512_set1_pd(4.6850901201934832188e+02), 
                                          _mm512_set1_pd(5.9169059852270512312e+03), 
                                          _mm512_set1_pd(2.4708152720399552679e+03)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[2]   = {_mm512_set1_pd(-2.4994418972832303646e+02), 
	                                  _mm512_set1_pd(2.1312714303849120380e+04)};
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d f[4]   = {_mm512_set1_pd(-1.6414452837299064100e+00),
	                                  _mm512_set1_pd(-2.9601657892958843866e+02), 
                                          _mm512_set1_pd(-1.7733784684952985886e+04),
                                          _mm512_set1_pd(-4.0320340761145482298e+05)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d g[3]   = {_mm512_set1_pd(-2.5064972445877992730e+02), 
	                                  _mm512_set1_pd(2.9865713163054025489e+04), 
                                          _mm512_set1_pd(-1.6128136304458193998e+06)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[10] = {_mm512_set1_pd(1.1394980557384778174e+02), 
	                                    _mm512_set1_pd(3.6832589957340267940e+03), 
                                            _mm512_set1_pd(3.1075408980684392399e+04), 
                                            _mm512_set1_pd(1.0577068948034021957e+05), 
                                            _mm512_set1_pd(1.7398867902565686251e+05), 
                                            _mm512_set1_pd(1.5097646353289914539e+05), 
                                            _mm512_set1_pd(7.1557062783764037541e+04), 
                                            _mm512_set1_pd(1.8321525870183537725e+04), 
                                            _mm512_set1_pd(2.3444738764199315021e+03), 
                                            _mm512_set1_pd(1.1600249425076035558e+02)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d qq[10] = {_mm512_set1_pd(2.0013443064949242491e+02), 
	                                    _mm512_set1_pd(4.4329628889746408858e+03), 
                                            _mm512_set1_pd(3.1474655750295278825e+04), 
                                            _mm512_set1_pd(9.7418829762268075784e+04), 
                                            _mm512_set1_pd(1.5144644673520157801e+05), 
                                            _mm512_set1_pd(1.2689839587977598727e+05), 
                                            _mm512_set1_pd(5.8824616785857027752e+04), 
                                            _mm512_set1_pd(1.4847228371802360957e+04), 
                                            _mm512_set1_pd(1.8821890840982713696e+03), 
                                            _mm512_set1_pd(9.2556599177304839811e+01)};

                          const __m512d one    = _mm512_set1_pd(1.0e+0);
                          const __m512d zero   = _mm512_set1_pd(0.0e+0);
                          const __m512d xsmall = _mm512_set1_pd(1.11e-16);
                          const __m512d xinf   = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax   = _mm512_set1_pd(705.342e+0);
                          register __m512d result,sumf,sumg,sump,sumq,temp;
                          register __m512d x,xx;
                          arg = _mm512_load_pd(&parg[0]);
                          x = arg;
                          if(_mm512_cmp_pd_mask(zero,x,_CMP_LT_OQ)) {
                             
                               if(_mm512_cmp_pd_mask(x,one,_CMP_LE_OQ)) {
                                   
                                   temp = xlog(x);
                                   
                                   if(_mm512_cmp_pd_mask(x,small,_CMP_LT_OQ)) {
                                       
                                       result = _mm512_sub_pd(_mm512_div_pd(p[5],p[1]),temp);
                                   }
                                   else {
                                       
                                        xx   = _mm512_mul_pd(x,x);
                                        sump = _mm512_fmadd_pd(pp[0],xx,
                                                      _mm512_fmadd_pd(pp[1],xx,
                                                            _mm512_fmadd_pd(pp[2],xx,
                                                                  _mm512_fmadd_pd(pp[3],xx,
                                                                        _mm512_fmadd_pd(pp[4],xx,p[5])))));
                                        sumq = _mm512_fmadd_pd(xx,q[0],
                                                      _mm512_add_pd(xx,q[1]));
                                        sumf = _mm512_fmadd_pd(f[0],xx,
                                                      _mm512_fmadd_pd(f[1],xx,
                                                             _mm512_fmadd_pd(f[2],xx,f[3])));
                                        sumg = _mm512_fmadd_pd(
                                                       _mm512_add_pd(xx,g[0]),xx,
                                                                 _mm512_fmadd_pd(g[1],xx,g[2]));
                                        register __m512d t0 = _mm512_div_pd(sump,sumq);
                                        register __m512d t1 = _mm512_sub_pd(_mm512_div_pd(temp,sumg),temp);
                                        result   = _mm512_sub_pd(t0,
                                                           _mm512_mul_pd(xx,
                                                                  _mm512_mul_pd(sumf,t1)));
                                        if(jint==2) {
                                            result = _mm512_mul_pd(result,xexp(x));
                                        }
                                   }
                               }
                               else if(jint==1 && 
                                       _mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                     
                                       result = zero;       
                               }
                               else {
                                       
                                      xx   = _mm512_div_pd(one,x);
                                      sumq = xx;
                                      sump = p[0];
                                      sump = _mm512_fmadd_pd(sump,xx,p[1]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[0]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[2]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[1]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[3]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[2]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[4]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[3]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[5]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[4]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[6]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[5]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[7]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[6]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[8]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[7]));
                                      sump = _mm512_fmadd_pd(sump,xx,p[9]);
                                      sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[8]));
                                      sumq = _mm512_add_pd(sumq,qq[9]);
                                      register __m512d t0 = _mm512_sqrt_pd(x);
                                      register __m512d t1 = _mm512_div_pd(sump,sumq);
                                      result = _mm512_div_pd(t1,t0);
                                      if(jint==1) {
                                         result = _mm512_mul_pd(result,xexp(negate_zmm8r8(x)));
                                      }
                               }
                          }  
                          else {
                                  result = xinf;
                          }  
                          
                          return (result);                  
	        }   
	        
	        
/*
!*****************************************************************************80
!
!! CALCK1 computes various K1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
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
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K1(x);
!    2, RESULT = exp(x) * K1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);
*/      


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calck1_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[5]  = {_mm512_set1_pd(4.8127070456878442310e-1), 
	                                   _mm512_set1_pd(9.9991373567429309922e+1), 
                                           _mm512_set1_pd(7.1885382604084798576e+3), 
                                           _mm512_set1_pd(1.7733324035147015630e+5), 
                                           _mm512_set1_pd(7.1938920065420586101e+5)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[3]  = {_mm512_set1_pd(-2.8143915754538725829e+2), 
	                                   _mm512_set1_pd(3.7264298672067697862e+4), 
                                           _mm512_set1_pd(-2.2149374878243304548e+6)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d f[5]  = {_mm512_set1_pd(-2.2795590826955002390e-1),
	                                   _mm512_set1_pd(-5.3103913335180275253e+1), 
                                           _mm512_set1_pd(-4.5051623763436087023e+3),
                                           _mm512_set1_pd(-1.4758069205414222471e+5), 
                                           _mm512_set1_pd(-1.3531161492785421328e+6)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d g[3]  = {_mm512_set1_pd(-3.0507151578787595807e+2), 
	                                   _mm512_set1_pd(4.3117653211351080007e+4), 
                                           _mm512_set1_pd(-2.7062322985570842656e+6)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[11] = {_mm512_set1_pd(6.4257745859173138767e-2),
	                                    _mm512_set1_pd(7.5584584631176030810e+0), 
                                            _mm512_set1_pd(1.3182609918569941308e+2), 
                                            _mm512_set1_pd(8.1094256146537402173e+2),
                                            _mm512_set1_pd(2.3123742209168871550e+3), 
                                            _mm512_set1_pd(3.4540675585544584407e+3), 
                                            _mm512_set1_pd(2.8590657697910288226e+3), 
                                            _mm512_set1_pd(1.3319486433183221990e+3), 
                                            _mm512_set1_pd(3.4122953486801312910e+2), 
                                            _mm512_set1_pd(4.4137176114230414036e+1), 
                                            _mm512_set1_pd(2.2196792496874548962e+0)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d qq[9] =  {_mm512_set1_pd(3.6001069306861518855e+1), 
	                                    _mm512_set1_pd(3.3031020088765390854e+2), 
                                            _mm512_set1_pd(1.2082692316002348638e+3), 
                                            _mm512_set1_pd(2.1181000487171943810e+3), 
                                            _mm512_set1_pd(1.9448440788918006154e+3), 
                                            _mm512_set1_pd(9.6929165726802648634e+2), 
                                            _mm512_set1_pd(2.5951223655579051357e+2), 
                                            _mm512_set1_pd(3.4552228452758912848e+1), 
                                            _mm512_set1_pd(1.7710478032601086579e+0)};
                          const __m512d one   = _mm512_set1_pd(1.0e+0);
                          const __m512d zero  = _mm512_set1_pd(0.0e+0);
                          const __m512d xleast= _mm512_set1_pd(2.23e-308);
                          const __m512d xsmall= _mm512_set1_pd(1.11e-16);
                          const __m512d xinf  = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax  = _mm512_set1_pd(705.343e+0);
                          register __m512d result,sumf,sumg,sump,sumq;
                          register __m512d x,xx,t0,t1,t2;
                          
                          x = arg;
                          if(_mm512_cmp_pd_mask(x,xleast,_CMP_LT_OQ)) {
                             
                              result = xinf;
                          }
                          else if(_mm512_cmp_pd_mask(x,one,_CMP_LE_OQ)) {
                              
                                if(_mm512_cmp_pd_mask(x,small,_CMP_LT_OQ)) {
                                   
                                    result = _mm512_div_pd(one,x);
                                }
                                else {
                                    
                                    xx   = _mm512_mul_pd(x,x);
                                    sump = _mm512_fmadd_pd(p[0],xx,
                                                      _mm512_fmadd_pd(p[1],xx,
                                                            _mm512_fmadd_pd(p[2],xx,
                                                                  _mm512_fmadd_pd(p[3],xx,
                                                                        _mm512_fmadd_pd(p[4],xx,q[2])))));
                                    sumq = _mm512_fmadd_pd(
                                                      _mm512_add_pd(xx,q[0]),xx,
                                                                _mm512_fmadd_pd(q[1],xx,q[2]));
                                    sumf = _mm512_fmadd_pd(f[0],xx,
                                                     _mm512_fmadd_pd(f[1],xx,
                                                              _mm512_fmadd_pd(f[2],xx,
                                                                       _mm512_fmadd_pd(f[3],xx,f[4]))));
                                    sumg = _mm512_fmadd_pd(
                                                      _mm512_add_pd(xx,g[0]),xx,
                                                                _mm512_fmadd_pd(g[1],xx,g[2]));
                                                                                  
                                    t0     = _mm512_div_pd(sumf,sumg);
                                    t2     = xlog(x);
                                    t1     = _mm512_div_pd(sump,sumq);
                                    result = _mm512_div_pd(_mm512_fmadd_pd(xx,t2,
                                                                    _mm512_add_pd(t0,t1)),x);
                                    if(jint==2) {
                                       result = _mm512_mul_pd(result,xexp(x));
                                    }                                             
                                }
                          } 
                          else if(jint==1 &&
                                  _mm512_cmp_pd_mask(xmamx,x,_CMP_LT_OQ)) {
                                   
                                     result = zero;
                          }    
                          else {
                          
                                 xx   = _mm512_div_pd(one,x);
                                 sumq = xx;
                                 sump = pp[0];
                                 sump = _mm512_fmadd_pd(sump,xx,p[1]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[0]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[2]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[1]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[3]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[2]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[4]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[3]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[5]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[4]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[6]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[5]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[7]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[6]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[8]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[7]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[9]);
                                 sump = _mm512_fmadd_pd(sump,xx,p[10]);
                                 sumq = _mm512_add_pd(sumq,qq[8]);
                                 t0   = _mm512_sqrt_pd(x);
                                 t1   = _mm512_div_pd(sumq,sumq);
                                 result= _mm512_div_pd(t1,t0);
                                 
                                 if(jint==1) {
                                    result = _mm512_mul_pd(result,
                                                  xexp(negate_zmm8r8(x)));
                                 }
                          } 
                          
                          return (result);                     
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calck1_zmm8r8_a(const double * __restrict __ATTR-ALIGN__(64) parg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d p[5]  = {_mm512_set1_pd(4.8127070456878442310e-1), 
	                                   _mm512_set1_pd(9.9991373567429309922e+1), 
                                           _mm512_set1_pd(7.1885382604084798576e+3), 
                                           _mm512_set1_pd(1.7733324035147015630e+5), 
                                           _mm512_set1_pd(7.1938920065420586101e+5)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[3]  = {_mm512_set1_pd(-2.8143915754538725829e+2), 
	                                   _mm512_set1_pd(3.7264298672067697862e+4), 
                                           _mm512_set1_pd(-2.2149374878243304548e+6)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d f[5]  = {_mm512_set1_pd(-2.2795590826955002390e-1),
	                                   _mm512_set1_pd(-5.3103913335180275253e+1), 
                                           _mm512_set1_pd(-4.5051623763436087023e+3),
                                           _mm512_set1_pd(-1.4758069205414222471e+5), 
                                           _mm512_set1_pd(-1.3531161492785421328e+6)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d g[3]  = {_mm512_set1_pd(-3.0507151578787595807e+2), 
	                                   _mm512_set1_pd(4.3117653211351080007e+4), 
                                           _mm512_set1_pd(-2.7062322985570842656e+6)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d pp[11] = {_mm512_set1_pd(6.4257745859173138767e-2),
	                                    _mm512_set1_pd(7.5584584631176030810e+0), 
                                            _mm512_set1_pd(1.3182609918569941308e+2), 
                                            _mm512_set1_pd(8.1094256146537402173e+2),
                                            _mm512_set1_pd(2.3123742209168871550e+3), 
                                            _mm512_set1_pd(3.4540675585544584407e+3), 
                                            _mm512_set1_pd(2.8590657697910288226e+3), 
                                            _mm512_set1_pd(1.3319486433183221990e+3), 
                                            _mm512_set1_pd(3.4122953486801312910e+2), 
                                            _mm512_set1_pd(4.4137176114230414036e+1), 
                                            _mm512_set1_pd(2.2196792496874548962e+0)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d qq[9] =  {_mm512_set1_pd(3.6001069306861518855e+1), 
	                                    _mm512_set1_pd(3.3031020088765390854e+2), 
                                            _mm512_set1_pd(1.2082692316002348638e+3), 
                                            _mm512_set1_pd(2.1181000487171943810e+3), 
                                            _mm512_set1_pd(1.9448440788918006154e+3), 
                                            _mm512_set1_pd(9.6929165726802648634e+2), 
                                            _mm512_set1_pd(2.5951223655579051357e+2), 
                                            _mm512_set1_pd(3.4552228452758912848e+1), 
                                            _mm512_set1_pd(1.7710478032601086579e+0)};
                          const __m512d one   = _mm512_set1_pd(1.0e+0);
                          const __m512d zero  = _mm512_set1_pd(0.0e+0);
                          const __m512d xleast= _mm512_set1_pd(2.23e-308);
                          const __m512d xsmall= _mm512_set1_pd(1.11e-16);
                          const __m512d xinf  = _mm512_set1_pd(1.79e+308);
                          const __m512d xmax  = _mm512_set1_pd(705.343e+0);
                          register __m512d result,sumf,sumg,sump,sumq;
                          register __m512d x,xx,t0,t1,t2;
                          
                          arg = _mm512_load_pd(&parg[0]);
                          x = arg;
                          if(_mm512_cmp_pd_mask(x,xleast,_CMP_LT_OQ)) {
                             
                              result = xinf;
                          }
                          else if(_mm512_cmp_pd_mask(x,one,_CMP_LE_OQ)) {
                              
                                if(_mm512_cmp_pd_mask(x,small,_CMP_LT_OQ)) {
                                   
                                    result = _mm512_div_pd(one,x);
                                }
                                else {
                                    
                                    xx   = _mm512_mul_pd(x,x);
                                    sump = _mm512_fmadd_pd(p[0],xx,
                                                      _mm512_fmadd_pd(p[1],xx,
                                                            _mm512_fmadd_pd(p[2],xx,
                                                                  _mm512_fmadd_pd(p[3],xx,
                                                                        _mm512_fmadd_pd(p[4],xx,q[2])))));
                                    sumq = _mm512_fmadd_pd(
                                                      _mm512_add_pd(xx,q[0]),xx,
                                                                _mm512_fmadd_pd(q[1],xx,q[2]));
                                    sumf = _mm512_fmadd_pd(f[0],xx,
                                                     _mm512_fmadd_pd(f[1],xx,
                                                              _mm512_fmadd_pd(f[2],xx,
                                                                       _mm512_fmadd_pd(f[3],xx,f[4]))));
                                    sumg = _mm512_fmadd_pd(
                                                      _mm512_add_pd(xx,g[0]),xx,
                                                                _mm512_fmadd_pd(g[1],xx,g[2]));
                                                                                  
                                    t0     = _mm512_div_pd(sumf,sumg);
                                    t2     = xlog(x);
                                    t1     = _mm512_div_pd(sump,sumq);
                                    result = _mm512_div_pd(_mm512_fmadd_pd(xx,t2,
                                                                    _mm512_add_pd(t0,t1)),x);
                                    if(jint==2) {
                                       result = _mm512_mul_pd(result,xexp(x));
                                    }                                             
                                }
                          } 
                          else if(jint==1 &&
                                  _mm512_cmp_pd_mask(xmamx,x,_CMP_LT_OQ)) {
                                   
                                     result = zero;
                          }    
                          else {
                          
                                 xx   = _mm512_div_pd(one,x);
                                 sumq = xx;
                                 sump = pp[0];
                                 sump = _mm512_fmadd_pd(sump,xx,p[1]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[0]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[2]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[1]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[3]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[2]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[4]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[3]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[5]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[4]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[6]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[5]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[7]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[6]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[8]);
                                 sumq = _mm512_mul_pd(_mm512_add_pd(sumq,qq[7]));
                                 sump = _mm512_fmadd_pd(sump,xx,p[9]);
                                 sump = _mm512_fmadd_pd(sump,xx,p[10]);
                                 sumq = _mm512_add_pd(sumq,qq[8]);
                                 t0   = _mm512_sqrt_pd(x);
                                 t1   = _mm512_div_pd(sumq,sumq);
                                 result= _mm512_div_pd(t1,t0);
                                 
                                 if(jint==1) {
                                    result = _mm512_mul_pd(result,
                                                  xexp(negate_zmm8r8(x)));
                                 }
                          } 
                          
                          return (result);                     
	         }
	         
	         
/*
  !*****************************************************************************80
!
!! CALERF computes various forms of the error function.
!
!  Discussion:
!
!    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
!    for a real argument x.
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
!  Reference:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Function,
!    Mathematics of Computation,
!    Volume 23, Number 107, July 1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT is 1, the
!    argument must be less than XBIG.  If JINT is 2, the argument
!    must lie between XNEG and XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = erf(x);
!    1, RESULT = erfc(x) = 1 - erf(x);
!    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, erf(x);
!    1, erfc(x);
!    2, exp(x*x)*erfc(x).
*/


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calerf_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d a[5]  = {_mm512_set1_pd(3.16112374387056560e+0),
	                                   _mm512_set1_pd(1.13864154151050156e+2), 
                                           _mm512_set1_pd(3.77485237685302021e+2),
                                           _mm512_set1_pd(3.20937758913846947e+3),
                                           _mm512_set1_pd(1.85777706184603153e-1)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d b[4]  = {_mm512_set1_pd(2.36012909523441209e+1),
	                                   _mm512_set1_pd(2.44024637934444173e+2), 
                                           _mm512_set1_pd(1.28261652607737228e+3),
                                           _mm512_set1_pd(2.84423683343917062e+3)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d c[9]  = {_mm512_set1_pd(5.64188496988670089e-1),
	                                   _mm512_set1_pd(8.88314979438837594e+0), 
                                           _mm512_set1_pd(6.61191906371416295e+1),
                                           _mm512_set1_pd(2.98635138197400131e+2), 
                                           _mm512_set1_pd(8.81952221241769090e+2),
                                           _mm512_set1_pd(1.71204761263407058e+3),
                                           _mm512_set1_pd(2.05107837782607147e+3),
                                           _mm512_set1_pd(1.23033935479799725e+3), 
                                           _mm512_set1_pd(2.15311535474403846e-8)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d d[8]  = {_mm512_set1_pd(1.57449261107098347e+1),
	                                   _mm512_set1_pd(1.17693950891312499e+2), 
                                           _mm512_set1_pd(5.37181101862009858e+2),
                                           _mm512_set1_pd(1.62138957456669019e+3), 
                                           _mm512_set1_pd(3.29079923573345963e+3),
                                           _mm512_set1_pd(4.36261909014324716e+3), 
                                           _mm512_set1_pd(3.43936767414372164e+3),
                                           _mm512_set1_pd(1.23033935480374942e+3)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d p[6]  = {_mm512_set1_pd(3.05326634961232344e-1),
	                                   _mm512_set1_pd(3.60344899949804439e-1), 
                                           _mm512_set1_pd(1.25781726111229246e-1),
                                           _mm512_set1_pd(1.60837851487422766e-2), 
                                           _mm512_set1_pd(6.58749161529837803e-4),
                                           _mm512_set1_pd(1.63153871373020978e-2)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[5]  = {_mm512_set1_pd(2.56852019228982242e+0),
	                                   _mm512_set1_pd(1.87295284992346047e+0), 
                                           _mm512_set1_pd(5.27905102951428412e-1),
                                           _mm512_set1_pd(6.05183413124413191e-2), 
                                           _mm512_set1_pd(2.33520497626869185e-3)};
                          const __m512d four =  _mm512_set1_pd(4.0e+0);
                          const __m512d one  =  _mm512_set1_pd(1.0e+0);
                          const __m512d half =  _mm512_set1_pd(0.5e+0);
                          const __m512d two  =  _mm512_set1_pd(2.0e+0);
                          const __m512d zero =  _mm512_set1_pd(0.0e+0);
                          const __m512d sqrpi=  _mm512_set1_pd(5.6418958354775628695e-1);
                          const __m512d thresh= _mm512_set1_pd(0.46875e+0);
                          const __m512d sixten= _mm512_set1_pd(16.0e+0);
                          const __m512d xinf  = _mm512_set1_pd(1.79e+308);
                          const __m512d xneg  = _mm512_set1_pd(-26.628e+0);
                          const __m512d xsmall= _mm512_set1_pd(1.11e-16);
                          const __m512d xbig  = _mm512_set1_pd(26.543e+0);
                          const __m512d xhuge = _mm512_set1_pd(6.71e+7);
                          const __m512d xmax  = _mm512_set1_pd(2.53e+307);   
                          register __m512d del,result,x,xden,xnum,y,ysq,t0,t1;
                          
                          x = arg;
                          y = _mm512_abs_pd(x);
                          
                          if(_mm512_cmp_pd_mask(y,thresh,_CMP_LE_OQ)) {
                             
                              ysq = zero;
                              if(_mm512_cmp_pd_mask(xsmall,y,_CMP_LT_OQ)) {
                                 ysq = _mm512_mul_pd(y,y);
                              }
                              
                              xnum = _mm512_mul_pd(a[4],ysq);
                              xden = ysq;
                              xnum = _mm512_mul_pd(_mm512_add_pd(xnum,a[0]),ysq);
                              xden = _mm512_mul_pd(_mm512_add_pd(xden,b[0]),ysq);
                              xnum = _mm512_mul_pd(_mm512_add_pd(xnum,a[1]),ysq);
                              xden = _mm512_mul_pd(_mm512_add_pd(xden,b[1]),ysq);
                              xnum = _mm512_mul_pd(_mm512_add_pd(xnum,a[2]),ysq);
                              xden = _mm512_mul_pd(_mm512_add_pd(xden,b[2]),ysq);
                              t0   = _mm512_add_pd(xnum,a[3]);
                              t1   = _mm512_add_pd(xden,b[3]);
                              result = _mm512_mul_pd(x,_mm512_div_pd(t0,t1));
                              
                              if(jint!=0) {
                                  result = _mm512_sub_pd(one,result);
                              }
                              if(jint==2) {
                                  result = _mm512_mul_pd(xexp(ysq),result);
                              }
                              return (result);
                          }  
                          else if(_mm512_cmp_pd(y,four,_CMP_LE_OQ)) {
                                  
                                  xnum = _mm512_mul_pd(c[8],y);
                                  xden = y;
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[0]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[0]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[1]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[1]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[2]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[2]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[3]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[3]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[4]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[4]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[5]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[5]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[6]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[6]),y);
                                  t0   = _mm512_add_pd(xnum,c[7]);
                                  t1   = _mm512_add_pd(xden,d[7]);
                                  result = _mm512_div_pd(t0,t1);
                                  
                                  if(jint==2) {
                                      __m512i ti = _mm512_cvt_roundpd_epi64(
                                                        _mm512_mul_pd(y,sixten),
                                                                         _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                                      ysq = _mm512_div_pd(_mm512_castsi512_pd(ti),sixten);
                                      del = _mm512_mul_pd(_mm512_sub_pd(y,ysq),
                                                          _mm512_add_pd(y,ysq));
                                      t0  = _mm512_mul_pd(negate_zmm8r8(ysq),ysq);
                                      t1  = xexp(t0);
                                      result = _mm512_mul_pd(t1,
                                             _mm512_mul_pd(negate_zmm8r8(del),result));
                                }
                          }   
                          else {
                          
                                 result = zero;
                                 
                                 if(_mm512_cmp_pd_mask(xbig,y,_CMP_LE_OQ)) {
                                      
                                      if(jint!=2 || 
                                          _mm512_cmp_pd_mask(xmax,y,_CMP_LE_OQ)) {
                                           goto L300;
                                      }
                                      
                                      if(_mm512_cmp_pd_mask(xhuge,y,_CMP_LE_OQ)) {
                                          result = _mm512_div_pd(sqrpi,y);
                                          goto L300;
                                      }
                                 }
                                 
                                 ysq = _mm512_div_pd(one,_mm512_mul_pd(y,y));
                                 xnum= _mm512_mul_pd(p[5],ysq);
                                 xden= ysq;
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[0]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[0]),ysq);
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[1]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[1]),ysq);
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[2]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[2]),ysq);
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[3]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[3]),ysq);
                                 t0  = _mm512_add_pd(xnum,p[4]);
                                 t1  = _mm512_add_pd(xden,q[4]);
                                 result = _mm512_mul_pd(ysq,_mm512_div_pd(t0,t1));
                                 result = _mm512_div_pd(_mm512_sub_pd(sqrpi,result),y);
                                 
                                 if(jint!=2) {
                                     __m512i ti = _mm512_cvt_roundpd_epi64(
                                                        _mm512_mul_pd(y,sixten),
                                                                         _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                                      ysq = _mm512_div_pd(_mm512_castsi512_pd(ti),sixten);
                                      del = _mm512_mul_pd(_mm512_sub_pd(y,ysq),
                                                          _mm512_add_pd(y,ysq));
                                      t0  = _mm512_mul_pd(negate_zmm8r8(ysq),ysq);
                                      t1  = xexp(t0);
                                      result = _mm512_mul_pd(t1,
                                             _mm512_mul_pd(negate_zmm8r8(del),result));
                                 }
                          } 
                          
                    L300:
                              
                              if(jint==0) {
                              
                                  result = _mm512_add_pd(_mm512_sub_pd(half,result),half);
                                  if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ)) {
                                     result = negate_zmm8r8(result);
                                  }
                                  
                              }
                              else if(jint==1) {
                                  
                                  if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ)) {
                                       result = _mm512_dub_pd(two,result);
                                  }
                                  
                              }
                              else {
                                 
                                  if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ)) {
                                     
                                      if(_mm512_cmp_pd_mask(x,xneg,_CMP_LT_OQ)) {
                                         
                                         result = xinf;
                                      }
                                      else {
                                        
                                           __m512i ti = _mm512_cvt_roundpd_epi64(
                                                        _mm512_mul_pd(y,sixten),
                                                                         _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                                           ysq = _mm512_div_pd(_mm512_castsi512_pd(ti),sixten);
                                           del = _mm512_mul_pd(_mm512_sub_pd(y,ysq),
                                                          _mm512_add_pd(y,ysq));
                                           y   = _mm512_mul_pd(_mm512_mul_pd(ysq,ysq),cexp(del));
                                           result = _mm512_sub_pd(_mm512_add_pd(y,y),result);
                                      }
                                  }
                              }
                              
                              return (result);
                                   
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d calerf_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                 const int32_t jint) {
	                                 
	                  __ATTR_ALIGN__(64) const static 
	                  __m512d a[5]  = {_mm512_set1_pd(3.16112374387056560e+0),
	                                   _mm512_set1_pd(1.13864154151050156e+2), 
                                           _mm512_set1_pd(3.77485237685302021e+2),
                                           _mm512_set1_pd(3.20937758913846947e+3),
                                           _mm512_set1_pd(1.85777706184603153e-1)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d b[4]  = {_mm512_set1_pd(2.36012909523441209e+1),
	                                   _mm512_set1_pd(2.44024637934444173e+2), 
                                           _mm512_set1_pd(1.28261652607737228e+3),
                                           _mm512_set1_pd(2.84423683343917062e+3)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d c[9]  = {_mm512_set1_pd(5.64188496988670089e-1),
	                                   _mm512_set1_pd(8.88314979438837594e+0), 
                                           _mm512_set1_pd(6.61191906371416295e+1),
                                           _mm512_set1_pd(2.98635138197400131e+2), 
                                           _mm512_set1_pd(8.81952221241769090e+2),
                                           _mm512_set1_pd(1.71204761263407058e+3),
                                           _mm512_set1_pd(2.05107837782607147e+3),
                                           _mm512_set1_pd(1.23033935479799725e+3), 
                                           _mm512_set1_pd(2.15311535474403846e-8)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d d[8]  = {_mm512_set1_pd(1.57449261107098347e+1),
	                                   _mm512_set1_pd(1.17693950891312499e+2), 
                                           _mm512_set1_pd(5.37181101862009858e+2),
                                           _mm512_set1_pd(1.62138957456669019e+3), 
                                           _mm512_set1_pd(3.29079923573345963e+3),
                                           _mm512_set1_pd(4.36261909014324716e+3), 
                                           _mm512_set1_pd(3.43936767414372164e+3),
                                           _mm512_set1_pd(1.23033935480374942e+3)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d p[6]  = {_mm512_set1_pd(3.05326634961232344e-1),
	                                   _mm512_set1_pd(3.60344899949804439e-1), 
                                           _mm512_set1_pd(1.25781726111229246e-1),
                                           _mm512_set1_pd(1.60837851487422766e-2), 
                                           _mm512_set1_pd(6.58749161529837803e-4),
                                           _mm512_set1_pd(1.63153871373020978e-2)};
                          __ATTR_ALIGN__(64) const static 
	                  __m512d q[5]  = {_mm512_set1_pd(2.56852019228982242e+0),
	                                   _mm512_set1_pd(1.87295284992346047e+0), 
                                           _mm512_set1_pd(5.27905102951428412e-1),
                                           _mm512_set1_pd(6.05183413124413191e-2), 
                                           _mm512_set1_pd(2.33520497626869185e-3)};
                          const __m512d four =  _mm512_set1_pd(4.0e+0);
                          const __m512d one  =  _mm512_set1_pd(1.0e+0);
                          const __m512d half =  _mm512_set1_pd(0.5e+0);
                          const __m512d two  =  _mm512_set1_pd(2.0e+0);
                          const __m512d zero =  _mm512_set1_pd(0.0e+0);
                          const __m512d sqrpi=  _mm512_set1_pd(5.6418958354775628695e-1);
                          const __m512d thresh= _mm512_set1_pd(0.46875e+0);
                          const __m512d sixten= _mm512_set1_pd(16.0e+0);
                          const __m512d xinf  = _mm512_set1_pd(1.79e+308);
                          const __m512d xneg  = _mm512_set1_pd(-26.628e+0);
                          const __m512d xsmall= _mm512_set1_pd(1.11e-16);
                          const __m512d xbig  = _mm512_set1_pd(26.543e+0);
                          const __m512d xhuge = _mm512_set1_pd(6.71e+7);
                          const __m512d xmax  = _mm512_set1_pd(2.53e+307);   
                          register __m512d del,result,x,xden,xnum,y,ysq,t0,t1;
                          
                          x = _mm512_load_pd(&parg[0]);
                          y = _mm512_abs_pd(x);
                          
                          if(_mm512_cmp_pd_mask(y,thresh,_CMP_LE_OQ)) {
                             
                              ysq = zero;
                              if(_mm512_cmp_pd_mask(xsmall,y,_CMP_LT_OQ)) {
                                 ysq = _mm512_mul_pd(y,y);
                              }
                              
                              xnum = _mm512_mul_pd(a[4],ysq);
                              xden = ysq;
                              xnum = _mm512_mul_pd(_mm512_add_pd(xnum,a[0]),ysq);
                              xden = _mm512_mul_pd(_mm512_add_pd(xden,b[0]),ysq);
                              xnum = _mm512_mul_pd(_mm512_add_pd(xnum,a[1]),ysq);
                              xden = _mm512_mul_pd(_mm512_add_pd(xden,b[1]),ysq);
                              xnum = _mm512_mul_pd(_mm512_add_pd(xnum,a[2]),ysq);
                              xden = _mm512_mul_pd(_mm512_add_pd(xden,b[2]),ysq);
                              t0   = _mm512_add_pd(xnum,a[3]);
                              t1   = _mm512_add_pd(xden,b[3]);
                              result = _mm512_mul_pd(x,_mm512_div_pd(t0,t1));
                              
                              if(jint!=0) {
                                  result = _mm512_sub_pd(one,result);
                              }
                              if(jint==2) {
                                  result = _mm512_mul_pd(xexp(ysq),result);
                              }
                              return (result);
                          }  
                          else if(_mm512_cmp_pd(y,four,_CMP_LE_OQ)) {
                                  
                                  xnum = _mm512_mul_pd(c[8],y);
                                  xden = y;
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[0]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[0]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[1]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[1]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[2]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[2]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[3]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[3]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[4]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[4]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[5]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[5]),y);
                                  xnum = _mm512_mul_pd(_mm512_add_pd(xnum,c[6]),y);
                                  xden = _mm512_mul_pd(_mm512_add_pd(xden,d[6]),y);
                                  t0   = _mm512_add_pd(xnum,c[7]);
                                  t1   = _mm512_add_pd(xden,d[7]);
                                  result = _mm512_div_pd(t0,t1);
                                  
                                  if(jint==2) {
                                      __m512i ti = _mm512_cvt_roundpd_epi64(
                                                        _mm512_mul_pd(y,sixten),
                                                                         _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                                      ysq = _mm512_div_pd(_mm512_castsi512_pd(ti),sixten);
                                      del = _mm512_mul_pd(_mm512_sub_pd(y,ysq),
                                                          _mm512_add_pd(y,ysq));
                                      t0  = _mm512_mul_pd(negate_zmm8r8(ysq),ysq);
                                      t1  = xexp(t0);
                                      result = _mm512_mul_pd(t1,
                                             _mm512_mul_pd(negate_zmm8r8(del),result));
                                }
                          }   
                          else {
                          
                                 result = zero;
                                 
                                 if(_mm512_cmp_pd_mask(xbig,y,_CMP_LE_OQ)) {
                                      
                                      if(jint!=2 || 
                                          _mm512_cmp_pd_mask(xmax,y,_CMP_LE_OQ)) {
                                           goto L300;
                                      }
                                      
                                      if(_mm512_cmp_pd_mask(xhuge,y,_CMP_LE_OQ)) {
                                          result = _mm512_div_pd(sqrpi,y);
                                          goto L300;
                                      }
                                 }
                                 
                                 ysq = _mm512_div_pd(one,_mm512_mul_pd(y,y));
                                 xnum= _mm512_mul_pd(p[5],ysq);
                                 xden= ysq;
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[0]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[0]),ysq);
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[1]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[1]),ysq);
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[2]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[2]),ysq);
                                 xnum= _mm512_mul_pd(_mm512_add_pd(xnum,p[3]),ysq);
                                 xden= _mm512_mul_pd(_mm512_add_pd(xden,q[3]),ysq);
                                 t0  = _mm512_add_pd(xnum,p[4]);
                                 t1  = _mm512_add_pd(xden,q[4]);
                                 result = _mm512_mul_pd(ysq,_mm512_div_pd(t0,t1));
                                 result = _mm512_div_pd(_mm512_sub_pd(sqrpi,result),y);
                                 
                                 if(jint!=2) {
                                     __m512i ti = _mm512_cvt_roundpd_epi64(
                                                        _mm512_mul_pd(y,sixten),
                                                                         _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                                      ysq = _mm512_div_pd(_mm512_castsi512_pd(ti),sixten);
                                      del = _mm512_mul_pd(_mm512_sub_pd(y,ysq),
                                                          _mm512_add_pd(y,ysq));
                                      t0  = _mm512_mul_pd(negate_zmm8r8(ysq),ysq);
                                      t1  = xexp(t0);
                                      result = _mm512_mul_pd(t1,
                                             _mm512_mul_pd(negate_zmm8r8(del),result));
                                 }
                          } 
                          
                    L300:
                              
                              if(jint==0) {
                              
                                  result = _mm512_add_pd(_mm512_sub_pd(half,result),half);
                                  if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ)) {
                                     result = negate_zmm8r8(result);
                                  }
                                  
                              }
                              else if(jint==1) {
                                  
                                  if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ)) {
                                       result = _mm512_dub_pd(two,result);
                                  }
                                  
                              }
                              else {
                                 
                                  if(_mm512_cmp_pd_mask(x,zero,_CMP_LT_OQ)) {
                                     
                                      if(_mm512_cmp_pd_mask(x,xneg,_CMP_LT_OQ)) {
                                         
                                         result = xinf;
                                      }
                                      else {
                                        
                                           __m512i ti = _mm512_cvt_roundpd_epi64(
                                                        _mm512_mul_pd(y,sixten),
                                                                         _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                                           ysq = _mm512_div_pd(_mm512_castsi512_pd(ti),sixten);
                                           del = _mm512_mul_pd(_mm512_sub_pd(y,ysq),
                                                          _mm512_add_pd(y,ysq));
                                           y   = _mm512_mul_pd(_mm512_mul_pd(ysq,ysq),cexp(del));
                                           result = _mm512_sub_pd(_mm512_add_pd(y,y),result);
                                      }
                                  }
                              }
                              
                              return (result);
                                   
	          }
	          
	          
	          
/*
    !*****************************************************************************80
!
!! CALJY0 computes various J0 and Y0 Bessel functions.
!
!  Discussion:
!
!    This routine computes zero-order Bessel functions of the first and
!    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!    for Y0, and |X| <= XMAX for J0.
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
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J0(x);
!    1, RESULT = Y0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, J0(x);
!    1, Y0(x);
*/


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d caljy0_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                 __ATTR_ALIGN__(64) const static 
	                 __m512d plg[4]  = {_mm512_set1_pd(-2.4562334077563243311e+1),
	                                    _mm512_set1_pd(2.3642701335621505212e+2), 
                                            _mm512_set1_pd(-5.4989956895857911039e+2),
                                            _mm512_set1_pd(3.5687548468071500413e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qlg[4]  = {_mm512_set1_pd(3.5553900764052419184e+1),
	                                    _mm512_set1_pd(1.9400230218539473193e+2), 
                                            _mm512_set1_pd(-3.3442903192607538956e+2),
                                            _mm512_set1_pd(1.7843774234035750207e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj0[7]  = {_mm512_set14_pd(6.6302997904833794242e+6),
	                                    _mm512_set1_pd(-6.2140700423540120665e+8), 
                                            _mm512_set1_pd(2.7282507878605942706e+10),
                                            _mm512_set1_pd(-4.1298668500990866786e+11), 
                                            _mm512_set1_pd(-1.2117036164593528341e-1), 
                                            _mm512_set1_pd(1.0344222815443188943e+2), 
                                            _mm512_set1_pd(-3.6629814655107086448e+4)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj0[5]  = {_mm512_set1_pd(4.5612696224219938200e+5), 
	                                    _mm512_set1_pd(1.3985097372263433271e+8), 
                                            _mm512_set1_pd(2.6328198300859648632e+10), 
                                            _mm512_set1_pd(2.3883787996332290397e+12), 
                                            _mm512_set1_pd(9.3614022392337710626e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj1[8]  = {_mm512_set1_pd(4.4176707025325087628e+3), 
	                                    _mm512_set1_pd(1.1725046279757103576e+4),
                                            _mm512_set1_pd(1.0341910641583726701e+4),
                                            _mm512_set1_pd(-7.2879702464464618998e+3), 
                                            _mm512_set1_pd(-1.2254078161378989535e+4),
                                            _mm512_set1_pd(-1.8319397969392084011e+3), 
                                            _mm512_set1_pd(4.8591703355916499363e+1), 
                                            _mm512_set1_pd(7.4321196680624245801e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj1[7]  = {_mm512_set1_pd(3.3307310774649071172e+2),
	                                    _mm512_set1_pd(-2.9458766545509337327e+3), 
                                            _mm512_set1_pd(1.8680990008359188352e+4),
                                            _mm512_set1_pd(-8.4055062591169562211e+4), 
                                            _mm512_set1_pd(2.4599102262586308984e+5),
                                            _mm512_set1_pd(-3.5783478026152301072e+5), 
                                            _mm512_set1_pd(-2.5258076240801555057e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py0[6]  = {_mm512_set1_pd(1.0102532948020907590e+4),
	                                    _mm512_set1_pd(-2.1287548474401797963e+6), 
                                            _mm512_set1_pd(2.0422274357376619816e+8),
                                            _mm512_set1_pd(-8.3716255451260504098e+9), 
                                            _mm512_set1_pd(1.0723538782003176831e+11),
                                            _mm512_set1_pd(-1.8402381979244993524e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy0[5]  = {_mm512_set1_pd(6.6475986689240190091e+2), 
	                                    _mm512_set1_pd(2.3889393209447253406e+5), 
                                            _mm512_set1_pd(5.5662956624278251596e+7), 
                                            _mm512_set1_pd(8.1617187777290363573e+9), 
                                            _mm512_set1_pd(5.8873865738997033405e+11)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py1[7]  = {_mm512_set1_pd(1.4566865832663635920e+4), 
	                                    _mm512_set1_pd(4.6905288611678631510e+6), 
                                            _mm512_set1_pd(-6.9590439394619619534e+8), 
                                            _mm512_set1_pd(4.3600098638603061642e+10), 
                                            _mm512_set1_pd(-5.5107435206722644429e+11),
                                            _mm512_set1_pd(-2.2213976967566192242e+13), 
                                            _mm512_set1_pd(1.7427031242901594547e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy1[6]  = {_mm512_set1_pd(8.3030857612070288823e+2), 
	                                    _mm512_set1_pd(4.0669982352539552018e+5), 
                                            _mm512_set1_pd(1.3960202770986831075e+8), 
                                            _mm512_set1_pd(3.4015103849971240096e+10), 
                                            _mm512_set1_pd(5.4266824419412347550e+12), 
                                            _mm512_set1_pd(4.3386146580707264428e+14)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py2[8]  = {_mm512_set1_pd(2.1363534169313901632e+4),
	                                    _mm512_set1_pd(-1.0085539923498211426e+7), 
                                            _mm512_set1_pd(2.1958827170518100757e+9),
                                            _mm512_set1_pd(-1.9363051266772083678e+11), 
                                            _mm512_set1_pd(-1.2829912364088687306e+11), 
                                            _mm512_set1_pd(6.7016641869173237784e+14), 
                                            _mm512_set1_pd(-8.0728726905150210443e+15),
                                            _mm512_set1_pd(-1.7439661319197499338e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy2[7]  = {_mm512_set1_pd(8.7903362168128450017e+2), 
	                                    _mm512_set1_pd(5.3924739209768057030e+5), 
                                            _mm512_set1_pd(2.4727219475672302327e+8), 
                                            _mm512_set1_pd(8.6926121104209825246e+10),
                                            _mm512_set1_pd(2.2598377924042897629e+13), 
                                            _mm512_set1_pd(3.9272425569640309819e+15),
                                            _mm512_set1_pd(3.4563724628846457519e+17)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p0[6]  =  {_mm512_set1_pd(3.4806486443249270347e+3), 
	                                    _mm512_set1_pd(2.1170523380864944322e+4), 
                                            _mm512_set1_pd(4.1345386639580765797e+4), 
                                            _mm512_set1_pd(2.2779090197304684302e+4),
                                            _mm512_set1_pd(8.8961548424210455236e-1), 
                                            _mm512_set1_pd(1.5376201909008354296e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q0[5]  =  {_mm512_set1_pd(3.5028735138235608207e+3), 
	                                    _mm512_set1_pd(2.1215350561880115730e+4), 
                                            _mm512_set1_pd(4.1370412495510416640e+4), 
                                            _mm512_set1_pd(2.2779090197304684318e+4),
                                            _mm512_set1_pd(1.5711159858080893649e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p1[6]  =  {_mm512_set1_pd(-2.2300261666214198472e+1),
	                                    _mm512_set1_pd(-1.1183429920482737611e+2), 
                                            _mm512_set1_pd(-1.8591953644342993800e+2),
                                            _mm512_set1_pd(-8.9226600200800094098e+1),
                                            _mm512_set1_pd(-8.8033303048680751817e+3),
                                            _mm512_set1_pd(-1.2441026745835638459e+00)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q1[5]  =  {_mm512_set1_pd(1.4887231232283756582e+3), 
	                                    _mm512_set1_pd(7.2642780169211018836e+3),
                                            _mm512_set1_pd(1.1951131543434613647e+4), 
                                            _mm512_set1_pd(5.7105024128512061905e+3),
                                            _mm512_set1_pd(9.0593769594993125859e+1)};
                         
                        const __m512d zero  = _mm512_set1_pd(0.0e+0);
                        const __m512d one   = _mm512_set1_pd(1.0e+0);
                        const __m512d three = _mm512_set1_pd(3.0e+0);
                        const __m512d four  = _mm512_set1_pd(4.0e+0);
                        const __m512d eight = _mm512_set1_pd(8.0e+0);
                        const __m512d five5 = _mm512_set1_pd(5.5e+0); 
                        const __m512d sixty4= _mm512_set1_pd(64.0e+0);
                        const __m512d oneov8= _mm512_set1_pd(0.125e+0); 
                        const __m512d p17   = _mm512_set1_pd(1.716e-1);
                        const __m512d two56 = _mm512_set1_pd(256.0e+0);
                        const __m512d cons  = _mm512_set1_pd(-1.1593151565841244881e-1);
                        const __m512d pi2   = _mm512_set1_pd(6.3661977236758134308e-1);
                        const __m512d twopi = _mm5512_set1_pd(6.2831853071795864769e+0);
                        const __m12d  twopi1= _mm12_set1_pd(6.28125e+0);
                        const __m512d twopi2= _mm512_set1_pd(1.9353071795864769253e-3);
                        const __m512d xmax  = _mm512_set1_pd(1.07e+09);
                        const __m512d xsmall= _mm512_set1_pd(9.31e-10);
                        const __m512d xinf  = _mm512_set1_pd(1.7e+38);
                        const __m512d xj0   = _mm512_set1_pd(2.4048255576957727686e+0);
                        const __m512d xj1   = _mm512_set1_pd(5.5200781102863106496e+0);
                        const __m512d xy0   = _mm512_set1_pd(8.9357696627916752158e-1);
                        const __m512d xy1   = _mm512_set1_pd(3.9576784193148578684e+0);
                        const __m512d xy2   = _mm512_set1_pd(7.0860510603017726976e+0);
                        const __m512d xj01  = _mm512_set1_pd(616.0e+0);
                        const __m512d xj02  = _mm512_set1_pd(-1.4244423042272313784e-3);
                        const __m512d xj11  = _mm512_set1_pd(1413.0e+0);
                        const __m512d xj12  = _mm512_set1_pd(5.4686028631064959660e-4);
                        const __m512d xy01  = _mm512_set1_pd(228.0e+0);
                        const __m512d xy02  = _mm512_set1_pd(2.9519662791675215849e-3);
                        const __m512d xy11  = _mm512_set1_pd(1013.0e+0);
                        const __m512d xy12  = _mm512_set1_pd(6.4716931485786837568e-4);
                        const __m512d xy21  = _mm512_set1_pd(1814.0e+0);
                        const __m512d xy22  = _mm512_set1_pd(1.1356030177269762362e-4);   
                        register __m512d ax,down,prod,resj;
                        register __m512d result,r0,r1,up,w;
                        register __m512d wsq,xden,xy,z,zsq;
                        register __m512d t0,t1;
                        
                        ax = _mm512_abs_pd(arg);
                        if(jint==1 &&
                                _mm512_cmp_pd_mask(arg,zero,_CMP_LE_OQ)) {
                              result = negate_zmm8r8(xinf);
                              return (result);       
                        }  
                        else if(_mm512_cmp_pd_mask(xmax,ax,_CMP_LT_OQ)) {
                              result = zero;
                              return (result);
                        }                   
                        if(_mm512_cmp_pd_mask(eight,ax,_CMP_LT_OQ)) {
                              goto L800;
                        }           
                        if(_mm512_cmp_pd_mask(ax,small,_CMP_LE_OQ)) {
                              if(jint==0) {
                                 result = one;
                              }
                              else {
                                 result = _mm512_fmadd_pd(pi2,xlog(ax),cons);
                                 return (result);
                              }
                        }
 /*
                            !  Calculate J0 for appropriate interval, preserving
                            !  accuracy near the zero of J0.
                        
 */                     
                        zsq = _mm512_mul_pd(ax,ax);
                        if(_mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ)) {
                            xnum = _mm512_fmadd_pd(pj0[4],zsq,
                                               _mm512_fmadd_pd(pj0[5],zsq,pj0[6]));
                            xden = _mm512_add_pd(zsq,qj0[4]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[0]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[3]);
                            t0   = _mm512_sub_pd(ax,_mm512_div_pd(xj01,two56));
                            t1   = _mm512_add_pd(ax,xj0);
                            prod = _mm512_mul_pd(_mm512_sub_pd(t0,xj02),t1);
                        }  
                        else {
                            wsq = _mm512_sub_pd(one,_mm512_div_pd(zsq,sixty4));
                            xnum= _mm512_fmadd_pd(pj1[6],wsq,pj1[7]);
                            xden= _mm512_add_pd(wsq,qj1[6]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[0]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[0]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[1]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[1]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[2]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[2]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[3]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[3]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[4]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[4]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[5]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[5]);
                            t0  = _mm512_sub_pd(ax,_mm512_div_pd(xj11,two56));
                            t1  = _mm512_add_pd(ax,xj1);
                            prod= _mm512_sub_pd(_mm512_mul_pd(t0,t1),xj12);
                       }
                       result = _mm512_mul_pd(prod,_mm512_div_pd(xnum,xden));
                       if(jint==0) {
                          return (result);
                       }
                       /*
                          Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
                          !  where xn is a zero of Y0.
                       */
                       
                       if(_mm512_cmp_pd_mask(ax,three,_CMP_LE_OQ)) {
                           up = _mm512_sub_pd(_mm512_sub_pd(ax,
                                              _mm512_div_pd(xy01,two56)),xy02);
                           xy = xy0;
                       }
                       else if(_mm512_cmp_pd_mask(ax,five5,_CMP_LE_OQ)) {
                           up = _mm512_sub_pd(_mm512_sub_pd(ax,
                                              _mm512_div_pd(xy11,two56)),xy12);
                           xy = xy1;
                       }
                       else {
                           up = _mm512_sub_pd(_mm512_sub_pd(ax,
                                              _mm512_div_pd(xy21,two56)),xy22);
                           xy = xy2; 
                       }
                       down = _mm512_mul_pd(ax,xy);
                       
                       if(_mm512_cmp_pd_mask(
                                       _mm512_abs_pd(up),
                                                 _mm512_mul_pd(p17,down),
                                                                  _CMP_LT_OQ)) {
                             w   = _mm512_div_pd(up,down);
                             wsq = _mm512_mul_pd(w,w);
                             xnum= plg[0];
                             xden= _mm512_add_pd(wsq,qlg[0]);
                             xnum= _mm512_fmadd_pd(xnum,wsq,plg[1]);
                             xden= _mm512_fmadd_pd(xden,wsq,qlg[1]);
                             xnum= _mm512_fmadd_pd(xnum,wsq,plg[2]);
                             xden= _mm512_fmadd_pd(xden,wsq,qlg[2]);
                             xnum= _mm512_fmadd_pd(xnum,wsq,plg[3]);
                             xden= _mm512_fmadd_pd(xden,wsq,qlg[3]); 
                             t0  = _mm512_div_pd(xnum,xden);
                             t1  = _mm512_mul_pd(pi,result);
                             resj= _mm512_mul_pd(t1,_mm512_mul_pd(w,t0));                                          
                       }
                       else {
                             t0  = _mm512_div_pd(xnum,xden);
                             t1  = _mm512_mul_pd(pi,result);
                             resj= _mm512_mul_pd(t1,xlog(t0));
                       }
                       
                       /*
                           Now calculate Y0 for appropriate interval, preserving
                           !  accuracy near the zero of Y0.
                       */
                       
                       if(_mm512_cmp_pd_mask(ax,three,_CMP_LE_OQ)) {
                           
                            xnum = _mm512_fmadd_pd(py0[5],zsq,py0[0]);
                            xden = _mm512_add_pd(zsq,qy0[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[3]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[4]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[4]);
                       }
                       else if(_mm512_cmp_pd_mask(ax,five5,_CMP_LE_OQ)) {
                            
                            xnum = _mm512_fmadd_pd(py1[6],zsq,py1[0]);
                            xden = _mm512_add_pd(zsq,qy1[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[3]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[4]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[4]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[5]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[5]);
                       }
                       else {
                         
                            xnum = _mm512_fmadd_pd(py2[7],zsq,py2[0]);
                            xden = _mm512_add_pd(zsq,qy2[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[3]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[4]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[4]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[5]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[5]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[6]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[6]);
                       }
                       t0 = _mm512_div_pd(xnum,xden);
                       t1 = _mm512_mul_pd(up,down);
                       result = _mm512_fmadd_pd(t0,t1,resj);
                       return (result);
               L800:
                       /*
                           Calculate J0 or Y0 for 8.0 < |ARG|.
                       */ 
                       z = _mm512_div_pd(eight,ax);
                       w = _mm512_div_pd(ax,twopi);
                       __m512i ti = mm512_cvt_roundpd_epi64(w,
                                                     _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                       w          = _mm512_add_pd(_mm512_castsi512_pd(ti),oneov8);
                       t0         = _mm512_mul_pd(w,twopi2);
                       t1         = _mm512_sub_pd(ax,w);
                       w          = _mm512_fmsub_pd(t1,twopi1,t0);
                       zsq        = _mm512_mul_pd(z,z);
                       xnum       = _mm512_fmadd_pd(p0[4],zsq,p0[5]);
                       xden       = _mm512_add_pd(zsq,q0[4]);
                       up         = _mm512_fmadd_pd(p1[4],zsq,p1[5]);
                       down       = _mm512_add_pd(zsq,q1[4]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[0]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[0]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[1]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[1]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[2]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[2]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[3]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[3]); 
                       r0         = _mm512_div_pd(xnum,xden);
                       r1         = _mm512_div_pd(up,down);
                       
                        t0 = _mm512_sqrt_pd(_mm512_div_pd(pi2,ax));
                        t1 = xcos(w);
                        register __m512d t2 = xsin(w);
                        register __m512d t3 = _mm512_mul_pd(z,r1);
                       if(jint==0) {
                          result = _mm512_fmsub_pd(t0,
                                               _mm512_mul_pd(r0,t1),
                                                             _mm512_mul_pd(t3,t2));
                       }  
                       else {
                          result = _mm512_fmadd_pd(t0,
                                               _mm512_mul_pd(r0,t2),
                                                             _mm512_mul_pd(t3,t1));
                       }  
                       
                       return (result);              
	         }
	         
	         
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d caljy0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                 const int32_t jint) {
	                                 
	                 __ATTR_ALIGN__(64) const static 
	                 __m512d plg[4]  = {_mm512_set1_pd(-2.4562334077563243311e+1),
	                                    _mm512_set1_pd(2.3642701335621505212e+2), 
                                            _mm512_set1_pd(-5.4989956895857911039e+2),
                                            _mm512_set1_pd(3.5687548468071500413e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qlg[4]  = {_mm512_set1_pd(3.5553900764052419184e+1),
	                                    _mm512_set1_pd(1.9400230218539473193e+2), 
                                            _mm512_set1_pd(-3.3442903192607538956e+2),
                                            _mm512_set1_pd(1.7843774234035750207e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj0[7]  = {_mm512_set14_pd(6.6302997904833794242e+6),
	                                    _mm512_set1_pd(-6.2140700423540120665e+8), 
                                            _mm512_set1_pd(2.7282507878605942706e+10),
                                            _mm512_set1_pd(-4.1298668500990866786e+11), 
                                            _mm512_set1_pd(-1.2117036164593528341e-1), 
                                            _mm512_set1_pd(1.0344222815443188943e+2), 
                                            _mm512_set1_pd(-3.6629814655107086448e+4)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj0[5]  = {_mm512_set1_pd(4.5612696224219938200e+5), 
	                                    _mm512_set1_pd(1.3985097372263433271e+8), 
                                            _mm512_set1_pd(2.6328198300859648632e+10), 
                                            _mm512_set1_pd(2.3883787996332290397e+12), 
                                            _mm512_set1_pd(9.3614022392337710626e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj1[8]  = {_mm512_set1_pd(4.4176707025325087628e+3), 
	                                    _mm512_set1_pd(1.1725046279757103576e+4),
                                            _mm512_set1_pd(1.0341910641583726701e+4),
                                            _mm512_set1_pd(-7.2879702464464618998e+3), 
                                            _mm512_set1_pd(-1.2254078161378989535e+4),
                                            _mm512_set1_pd(-1.8319397969392084011e+3), 
                                            _mm512_set1_pd(4.8591703355916499363e+1), 
                                            _mm512_set1_pd(7.4321196680624245801e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj1[7]  = {_mm512_set1_pd(3.3307310774649071172e+2),
	                                    _mm512_set1_pd(-2.9458766545509337327e+3), 
                                            _mm512_set1_pd(1.8680990008359188352e+4),
                                            _mm512_set1_pd(-8.4055062591169562211e+4), 
                                            _mm512_set1_pd(2.4599102262586308984e+5),
                                            _mm512_set1_pd(-3.5783478026152301072e+5), 
                                            _mm512_set1_pd(-2.5258076240801555057e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py0[6]  = {_mm512_set1_pd(1.0102532948020907590e+4),
	                                    _mm512_set1_pd(-2.1287548474401797963e+6), 
                                            _mm512_set1_pd(2.0422274357376619816e+8),
                                            _mm512_set1_pd(-8.3716255451260504098e+9), 
                                            _mm512_set1_pd(1.0723538782003176831e+11),
                                            _mm512_set1_pd(-1.8402381979244993524e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy0[5]  = {_mm512_set1_pd(6.6475986689240190091e+2), 
	                                    _mm512_set1_pd(2.3889393209447253406e+5), 
                                            _mm512_set1_pd(5.5662956624278251596e+7), 
                                            _mm512_set1_pd(8.1617187777290363573e+9), 
                                            _mm512_set1_pd(5.8873865738997033405e+11)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py1[7]  = {_mm512_set1_pd(1.4566865832663635920e+4), 
	                                    _mm512_set1_pd(4.6905288611678631510e+6), 
                                            _mm512_set1_pd(-6.9590439394619619534e+8), 
                                            _mm512_set1_pd(4.3600098638603061642e+10), 
                                            _mm512_set1_pd(-5.5107435206722644429e+11),
                                            _mm512_set1_pd(-2.2213976967566192242e+13), 
                                            _mm512_set1_pd(1.7427031242901594547e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy1[6]  = {_mm512_set1_pd(8.3030857612070288823e+2), 
	                                    _mm512_set1_pd(4.0669982352539552018e+5), 
                                            _mm512_set1_pd(1.3960202770986831075e+8), 
                                            _mm512_set1_pd(3.4015103849971240096e+10), 
                                            _mm512_set1_pd(5.4266824419412347550e+12), 
                                            _mm512_set1_pd(4.3386146580707264428e+14)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py2[8]  = {_mm512_set1_pd(2.1363534169313901632e+4),
	                                    _mm512_set1_pd(-1.0085539923498211426e+7), 
                                            _mm512_set1_pd(2.1958827170518100757e+9),
                                            _mm512_set1_pd(-1.9363051266772083678e+11), 
                                            _mm512_set1_pd(-1.2829912364088687306e+11), 
                                            _mm512_set1_pd(6.7016641869173237784e+14), 
                                            _mm512_set1_pd(-8.0728726905150210443e+15),
                                            _mm512_set1_pd(-1.7439661319197499338e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy2[7]  = {_mm512_set1_pd(8.7903362168128450017e+2), 
	                                    _mm512_set1_pd(5.3924739209768057030e+5), 
                                            _mm512_set1_pd(2.4727219475672302327e+8), 
                                            _mm512_set1_pd(8.6926121104209825246e+10),
                                            _mm512_set1_pd(2.2598377924042897629e+13), 
                                            _mm512_set1_pd(3.9272425569640309819e+15),
                                            _mm512_set1_pd(3.4563724628846457519e+17)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p0[6]  =  {_mm512_set1_pd(3.4806486443249270347e+3), 
	                                    _mm512_set1_pd(2.1170523380864944322e+4), 
                                            _mm512_set1_pd(4.1345386639580765797e+4), 
                                            _mm512_set1_pd(2.2779090197304684302e+4),
                                            _mm512_set1_pd(8.8961548424210455236e-1), 
                                            _mm512_set1_pd(1.5376201909008354296e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q0[5]  =  {_mm512_set1_pd(3.5028735138235608207e+3), 
	                                    _mm512_set1_pd(2.1215350561880115730e+4), 
                                            _mm512_set1_pd(4.1370412495510416640e+4), 
                                            _mm512_set1_pd(2.2779090197304684318e+4),
                                            _mm512_set1_pd(1.5711159858080893649e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p1[6]  =  {_mm512_set1_pd(-2.2300261666214198472e+1),
	                                    _mm512_set1_pd(-1.1183429920482737611e+2), 
                                            _mm512_set1_pd(-1.8591953644342993800e+2),
                                            _mm512_set1_pd(-8.9226600200800094098e+1),
                                            _mm512_set1_pd(-8.8033303048680751817e+3),
                                            _mm512_set1_pd(-1.2441026745835638459e+00)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q1[5]  =  {_mm512_set1_pd(1.4887231232283756582e+3), 
	                                    _mm512_set1_pd(7.2642780169211018836e+3),
                                            _mm512_set1_pd(1.1951131543434613647e+4), 
                                            _mm512_set1_pd(5.7105024128512061905e+3),
                                            _mm512_set1_pd(9.0593769594993125859e+1)};
                         
                        const __m512d zero  = _mm512_set1_pd(0.0e+0);
                        const __m512d one   = _mm512_set1_pd(1.0e+0);
                        const __m512d three = _mm512_set1_pd(3.0e+0);
                        const __m512d four  = _mm512_set1_pd(4.0e+0);
                        const __m512d eight = _mm512_set1_pd(8.0e+0);
                        const __m512d five5 = _mm512_set1_pd(5.5e+0); 
                        const __m512d sixty4= _mm512_set1_pd(64.0e+0);
                        const __m512d oneov8= _mm512_set1_pd(0.125e+0); 
                        const __m512d p17   = _mm512_set1_pd(1.716e-1);
                        const __m512d two56 = _mm512_set1_pd(256.0e+0);
                        const __m512d cons  = _mm512_set1_pd(-1.1593151565841244881e-1);
                        const __m512d pi2   = _mm512_set1_pd(6.3661977236758134308e-1);
                        const __m512d twopi = _mm5512_set1_pd(6.2831853071795864769e+0);
                        const __m12d  twopi1= _mm12_set1_pd(6.28125e+0);
                        const __m512d twopi2= _mm512_set1_pd(1.9353071795864769253e-3);
                        const __m512d xmax  = _mm512_set1_pd(1.07e+09);
                        const __m512d xsmall= _mm512_set1_pd(9.31e-10);
                        const __m512d xinf  = _mm512_set1_pd(1.7e+38);
                        const __m512d xj0   = _mm512_set1_pd(2.4048255576957727686e+0);
                        const __m512d xj1   = _mm512_set1_pd(5.5200781102863106496e+0);
                        const __m512d xy0   = _mm512_set1_pd(8.9357696627916752158e-1);
                        const __m512d xy1   = _mm512_set1_pd(3.9576784193148578684e+0);
                        const __m512d xy2   = _mm512_set1_pd(7.0860510603017726976e+0);
                        const __m512d xj01  = _mm512_set1_pd(616.0e+0);
                        const __m512d xj02  = _mm512_set1_pd(-1.4244423042272313784e-3);
                        const __m512d xj11  = _mm512_set1_pd(1413.0e+0);
                        const __m512d xj12  = _mm512_set1_pd(5.4686028631064959660e-4);
                        const __m512d xy01  = _mm512_set1_pd(228.0e+0);
                        const __m512d xy02  = _mm512_set1_pd(2.9519662791675215849e-3);
                        const __m512d xy11  = _mm512_set1_pd(1013.0e+0);
                        const __m512d xy12  = _mm512_set1_pd(6.4716931485786837568e-4);
                        const __m512d xy21  = _mm512_set1_pd(1814.0e+0);
                        const __m512d xy22  = _mm512_set1_pd(1.1356030177269762362e-4);   
                        register __m512d ax,down,prod,resj;
                        register __m512d result,r0,r1,up,w;
                        register __m512d wsq,xden,xy,z,zsq;
                        register __m512d t0,t1;
                        arg = _mm512_load_pd(&parg[0]);
                        ax  = _mm512_abs_pd(arg);
                        if(jint==1 &&
                                _mm512_cmp_pd_mask(arg,zero,_CMP_LE_OQ)) {
                              result = negate_zmm8r8(xinf);
                              return (result);       
                        }  
                        else if(_mm512_cmp_pd_mask(xmax,ax,_CMP_LT_OQ)) {
                              result = zero;
                              return (result);
                        }                   
                        if(_mm512_cmp_pd_mask(eight,ax,_CMP_LT_OQ)) {
                              goto L800;
                        }           
                        if(_mm512_cmp_pd_mask(ax,small,_CMP_LE_OQ)) {
                              if(jint==0) {
                                 result = one;
                              }
                              else {
                                 result = _mm512_fmadd_pd(pi2,xlog(ax),cons);
                                 return (result);
                              }
                        }
 /*
                            !  Calculate J0 for appropriate interval, preserving
                            !  accuracy near the zero of J0.
                        
 */                     
                        zsq = _mm512_mul_pd(ax,ax);
                        if(_mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ)) {
                            xnum = _mm512_fmadd_pd(pj0[4],zsq,
                                               _mm512_fmadd_pd(pj0[5],zsq,pj0[6]));
                            xden = _mm512_add_pd(zsq,qj0[4]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[0]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,pj0[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qj0[3]);
                            t0   = _mm512_sub_pd(ax,_mm512_div_pd(xj01,two56));
                            t1   = _mm512_add_pd(ax,xj0);
                            prod = _mm512_mul_pd(_mm512_sub_pd(t0,xj02),t1);
                        }  
                        else {
                            wsq = _mm512_sub_pd(one,_mm512_div_pd(zsq,sixty4));
                            xnum= _mm512_fmadd_pd(pj1[6],wsq,pj1[7]);
                            xden= _mm512_add_pd(wsq,qj1[6]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[0]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[0]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[1]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[1]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[2]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[2]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[3]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[3]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[4]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[4]);
                            xnum= _mm512_fmadd_pd(xnum,wsq,pj1[5]);
                            xden= _mm512_fmadd_pd(xden,wsq,qj1[5]);
                            t0  = _mm512_sub_pd(ax,_mm512_div_pd(xj11,two56));
                            t1  = _mm512_add_pd(ax,xj1);
                            prod= _mm512_sub_pd(_mm512_mul_pd(t0,t1),xj12);
                       }
                       result = _mm512_mul_pd(prod,_mm512_div_pd(xnum,xden));
                       if(jint==0) {
                          return (result);
                       }
                       /*
                          Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
                          !  where xn is a zero of Y0.
                       */
                       
                       if(_mm512_cmp_pd_mask(ax,three,_CMP_LE_OQ)) {
                           up = _mm512_sub_pd(_mm512_sub_pd(ax,
                                              _mm512_div_pd(xy01,two56)),xy02);
                           xy = xy0;
                       }
                       else if(_mm512_cmp_pd_mask(ax,five5,_CMP_LE_OQ)) {
                           up = _mm512_sub_pd(_mm512_sub_pd(ax,
                                              _mm512_div_pd(xy11,two56)),xy12);
                           xy = xy1;
                       }
                       else {
                           up = _mm512_sub_pd(_mm512_sub_pd(ax,
                                              _mm512_div_pd(xy21,two56)),xy22);
                           xy = xy2; 
                       }
                       down = _mm512_mul_pd(ax,xy);
                       
                       if(_mm512_cmp_pd_mask(
                                       _mm512_abs_pd(up),
                                                 _mm512_mul_pd(p17,down),
                                                                  _CMP_LT_OQ)) {
                             w   = _mm512_div_pd(up,down);
                             wsq = _mm512_mul_pd(w,w);
                             xnum= plg[0];
                             xden= _mm512_add_pd(wsq,qlg[0]);
                             xnum= _mm512_fmadd_pd(xnum,wsq,plg[1]);
                             xden= _mm512_fmadd_pd(xden,wsq,qlg[1]);
                             xnum= _mm512_fmadd_pd(xnum,wsq,plg[2]);
                             xden= _mm512_fmadd_pd(xden,wsq,qlg[2]);
                             xnum= _mm512_fmadd_pd(xnum,wsq,plg[3]);
                             xden= _mm512_fmadd_pd(xden,wsq,qlg[3]); 
                             t0  = _mm512_div_pd(xnum,xden);
                             t1  = _mm512_mul_pd(pi,result);
                             resj= _mm512_mul_pd(t1,_mm512_mul_pd(w,t0));                                          
                       }
                       else {
                             t0  = _mm512_div_pd(xnum,xden);
                             t1  = _mm512_mul_pd(pi,result);
                             resj= _mm512_mul_pd(t1,xlog(t0));
                       }
                       
                       /*
                           Now calculate Y0 for appropriate interval, preserving
                           !  accuracy near the zero of Y0.
                       */
                       
                       if(_mm512_cmp_pd_mask(ax,three,_CMP_LE_OQ)) {
                           
                            xnum = _mm512_fmadd_pd(py0[5],zsq,py0[0]);
                            xden = _mm512_add_pd(zsq,qy0[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[3]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py0[4]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy0[4]);
                       }
                       else if(_mm512_cmp_pd_mask(ax,five5,_CMP_LE_OQ)) {
                            
                            xnum = _mm512_fmadd_pd(py1[6],zsq,py1[0]);
                            xden = _mm512_add_pd(zsq,qy1[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[3]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[4]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[4]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py1[5]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy1[5]);
                       }
                       else {
                         
                            xnum = _mm512_fmadd_pd(py2[7],zsq,py2[0]);
                            xden = _mm512_add_pd(zsq,qy2[0]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[1]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[1]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[2]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[2]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[3]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[3]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[4]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[4]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[5]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[5]);
                            xnum = _mm512_fmadd_pd(xnum,zsq,py2[6]);
                            xden = _mm512_fmadd_pd(xden,zsq,qy2[6]);
                       }
                       t0 = _mm512_div_pd(xnum,xden);
                       t1 = _mm512_mul_pd(up,down);
                       result = _mm512_fmadd_pd(t0,t1,resj);
                       return (result);
               L800:
                       /*
                           Calculate J0 or Y0 for 8.0 < |ARG|.
                       */ 
                       z = _mm512_div_pd(eight,ax);
                       w = _mm512_div_pd(ax,twopi);
                       __m512i ti = mm512_cvt_roundpd_epi64(w,
                                                     _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                       w          = _mm512_add_pd(_mm512_castsi512_pd(ti),oneov8);
                       t0         = _mm512_mul_pd(w,twopi2);
                       t1         = _mm512_sub_pd(ax,w);
                       w          = _mm512_fmsub_pd(t1,twopi1,t0);
                       zsq        = _mm512_mul_pd(z,z);
                       xnum       = _mm512_fmadd_pd(p0[4],zsq,p0[5]);
                       xden       = _mm512_add_pd(zsq,q0[4]);
                       up         = _mm512_fmadd_pd(p1[4],zsq,p1[5]);
                       down       = _mm512_add_pd(zsq,q1[4]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[0]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[0]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[1]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[1]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[2]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[2]);
                       xnum       = _mm512_fmadd_pd(xnum,zsq,p0[3]);
                       xden       = _mm512_fmadd_pd(xden,zsq,q0[3]); 
                       r0         = _mm512_div_pd(xnum,xden);
                       r1         = _mm512_div_pd(up,down);
                       
                        t0 = _mm512_sqrt_pd(_mm512_div_pd(pi2,ax));
                        t1 = xcos(w);
                        register __m512d t2 = xsin(w);
                        register __m512d t3 = _mm512_mul_pd(z,r1);
                       if(jint==0) {
                          result = _mm512_fmsub_pd(t0,
                                               _mm512_mul_pd(r0,t1),
                                                             _mm512_mul_pd(t3,t2));
                       }  
                       else {
                          result = _mm512_fmadd_pd(t0,
                                               _mm512_mul_pd(r0,t2),
                                                             _mm512_mul_pd(t3,t1));
                       }  
                       
                       return (result);              
	         }
	         
	         
	         /*
	             !*****************************************************************************80
!
!! CALJY1 computes various J1 and Y1 Bessel functions.
!
!  Discussion:
!
!    This routine computes first-order Bessel functions of the first and
!    second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!    for Y1, and |X| <= XMAX for J1.
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
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J1(x);
!    1, RESULT = Y1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, J1(x);
!    1, Y1(x);
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d caljy1_zmm8r8(const __m512d arg,
	                                 const int32_t jint) {
	                                 
	                 __ATTR_ALIGN__(64) const static 
	                 __m512d plg[4]  = {_mm512_set1_pd(-2.4562334077563243311e+1),
	                                    _mm512_set1_pd(2.3642701335621505212e+2),
                                            _mm512_set1_pd(-5.4989956895857911039e+2),
                                            _mm512_set1_pd(3.5687548468071500413e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qlg[4]  = {_mm512_set1_pd(-3.5553900764052419184e+1),
	                                    _mm512_set1_pd(1.9400230218539473193e+2),
                                            _mm512_set1_pd(-3.3442903192607538956e+2),
                                            _mm512_set1_pd(1.7843774234035750207e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj0[7]  = {_mm512_set1_pd(9.8062904098958257677e+5,
	                                    _mm512_set1_pd(-1.1548696764841276794e+8), 
                                            _mm512_set1_pd(6.6781041261492395835e+9),
                                            _mm512_set1_pd(-1.4258509801366645672e+11), 
                                            _mm512_set1_pd(-4.4615792982775076130e+3), 
                                            _mm512_set1_pd(1.0650724020080236441e+1),
                                            _mm512_set1_pd(-1.0767857011487300348e-2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj0[5]  = {_mm512_set1_pd(5.9117614494174794095e+5), 
	                                    _mm512_set1_pd(2.0228375140097033958e+8), 
                                            _mm512_set1_pd(4.2091902282580133541e+10), 
                                            _mm512_set1_pd(4.1868604460820175290e+12), 
                                            _mm512_set1_pd(1.0742272239517380498e+03)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj1[8]  = {_mm512_set1_pd(4.6179191852758252280e+00),
	                                    _mm512_set1_pd(-7.1329006872560947377e+3),
                                            _mm512_set1_pd(4.5039658105749078904e+6),
                                            _mm512_set1_pd(-1.4437717718363239107e+9),
                                            _mm512_set1_pd(2.3569285397217157313e+11),
                                            _mm512_set1_pd(-1.6324168293282543629e+13),
                                            _mm512_set1_pd(1.1357022719979468624e+14), 
                                            _mm512_set1_pd(1.0051899717115285432e+15)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj1[7]  = {_mm512_set1_pd(1.1267125065029138050e+6), 
	                                    _mm512_set1_pd(6.4872502899596389593e+8),
                                            _mm512_set1_pd(2.7622777286244082666e+11), 
                                            _mm512_set1_pd(8.4899346165481429307e+13),
                                            _mm512_set1_pd(1.7128800897135812012e+16), 
                                            _mm512_set1_pd(1.7253905888447681194e+18), 
                                            _mm512_set1_pd(1.3886978985861357615e+3)}; 
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py0[7]  = {_mm512_set1_pd(2.2157953222280260820e+5),
	                                    _mm512_set1_pd(-5.9157479997408395984e+7), 
                                            _mm512_set1_pd(7.2144548214502560419e+9),
                                            _mm512_set1_pd(-3.7595974497819597599e+11),
                                            _mm512_set1_pd(5.4708611716525426053e+12), 
                                            _mm512_set1_pd(4.0535726612579544093e+13), 
                                            _mm512_set1_pd(-3.1714424660046133456e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy0[6]  = {_mm512_set1_pd(8.2079908168393867438e+2), 
	                                    _mm512_set1_pd(3.8136470753052572164e+5),
                                            _mm512_set1_pd(1.2250435122182963220e+8), 
                                            _mm512_set1_pd(2.7800352738690585613e+10),
                                            _mm512_set1_pd(4.1272286200406461981e+12), 
                                            _mm512_set_pd(3.0737873921079286084e+14)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py1[9]  =  {_mm512_set1_pd(1.9153806858264202986e+6),
	                                     _mm512_set1_pd(-1.1957961912070617006e+9), 
                                             _mm512_set1_pd(3.7453673962438488783e+11),
                                             _mm512_set1_pd(-5.9530713129741981618e+13), 
                                             _mm512_set1_pd(4.0686275289804744814e+15),
                                             _mm512_set1_pd(-2.3638408497043134724e+16),
                                             _mm512_set1_pd(-5.6808094574724204577e+18), 
                                             _mm512_set1_pd(1.1514276357909013326e+19), 
                                             _mm512_set1_pd(-1.2337180442012953128e+3)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy1[8]  =  {_mm512_set1_pd(1.2855164849321609336e+3), 
	                                     _mm512_set1_pd(1.0453748201934079734e+6), 
                                             _mm512_set1_pd(6.3550318087088919566e+8), 
                                             _mm512_set1_pd(3.0221766852960403645e+11), 
                                             _mm512_set1_pd(1.1187010065856971027e+14), 
                                             _mm512_set1_pd(3.0837179548112881950e+16),
                                             _mm512_set1_pd(5.6968198822857178911e+18), 
                                             _mm512_set1_pd(5.3321844313316185697e+20)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p0[6]  =   {_mm512_set1_pd(-1.0982405543459346727e+5),
	                                     _mm512_set1_pd(-1.5235293511811373833e+6),
                                             _mm512_set1_pd(-6.6033732483649391093e06),
                                             _mm512_set1_pd(-9.9422465050776411957e+6),
                                             _mm512_set1_pd(-4.4357578167941278571e+6),
                                             _mm512_set1_pd(-1.6116166443246101165e+3)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q0[6]  =   {_mm512_set1_pd(-1.0726385991103820119e+5),
	                                     _mm512_set1_pd(-1.5118095066341608816e+6),
                                             _mm512_set1_pd(-6.5853394797230870728e+6),
                                             _mm512_set1_pd(-9.9341243899345856590e+6), 
                                             _mm512_set1_pd(-4.4357578167941278568e+6),
                                             _mm512_set1_pd(-1.4550094401904961825e+3)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p1[6]  =   {_mm512_set1_pd(1.7063754290207680021e+3), 
	                                     _mm512_set1_pd(1.8494262873223866797e+4), 
                                             _mm512_set1_pd(6.6178836581270835179e+4), 
                                             _mm512_set1_pd(8.5145160675335701966e+4),
                                             _mm512_set1_pd(3.3220913409857223519e+4), 
                                             _mm512_set1_pd(3.5265133846636032186e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q1[6]  =   {_mm512_set1_pd(3.7890229745772202641e+4), 
	                                     _mm512_set1_pd(4.0029443582266975117e+5),
                                             _mm512_set1_pd(1.4194606696037208929e+6), 
                                             _mm512_set1_pd(1.8194580422439972989e+6),
                                             _mm512_set1_pd(7.0871281941028743574e+5), 
                                             _mm512_set1_pd(8.6383677696049909675e+2)};
                         const __m512d eight = _mm512_set1_pd(8.0e+0);
                         const __m512d four  = _mm512_set1_pd(4.0e+0);
                         const __m512d  half = _mm512_set1_pd(0.5);
                         const __m512d throv8= _mm512_set1_pd(0.375);
                         const __m512d pi2   = _mm512_set1_pd(6.3661977236758134308e-1);
                         const __m512d p17   = _mm512_set1_pd(1.716e-1);
                         const __m512d twopi = _mm512_set1_pd(6.2831853071795864769e+0);
                         const __m512d zero  = _mm512_set1_pd(0.0);
                         const __m512d twopi1= _mm512_set1_pd(6.28125e+0);
                         const __m512d twopi2= _mm512_set1_pd(1.9353071795864769253e-3);
                         const __m512d two56 = _mm512_set1_pd(256.0e+0);
                         const __m512d rtpi2 = _mm512_set1_pd(7.9788456080286535588e-1);
                         const __m512d xmax  = _mm512_set1_pd(1.07e+9);
                         const __m512d xsmall= _mm512_set1_pd(9.31e-10);
                         const __m512d xinf  = _mm512_set1_pd(1.7e+38);  
                         const __m512d xj0   = _mm12_set1_pd(3.8317059702075123156e+0);
                         const __m512d xj1   = _mm512_set1_pd(7.0155866698156187535e+0);
                         const __m512d xy0   = _mm512_set1_pd(2.1971413260310170351e+0);
                         const __m512d xy1   = _mm512_set1_pd(5.4296810407941351328e+0);
                         const __m512d xj01  = _mm512_set1_pd(981.0e+0);
                         const __m512d xj02  = _mm512_set1_pd(-3.2527979248768438556e-4);
                         const __m512d xj11  = _mm512_set1_pd(1796.0e+0);
                         const __m512d xj12  = _mm512_set1_pd(-3.8330184381246462950e-5);
                         const __m512d xy01  = _mm512_set1_pd(562.0e+0);
                         const __m512d xy02  = _mm512_set1_pd(1.8288260310170351490e-3);
                         const __m512d xy11  = _mm512_set1_pd(1390.0e+0);
                         const __m512d xy12  = _mm512_set1_pd(-6.4592058648672279948e-6);
                         register __m512d ax,down,prod,resj,result;
                         register __m512d r0,r1,up,w,wsq;
                         register __m512d xden,xnum,t0,t1,z,zsq;
                         
                         ax = _mm512_abs_pd(arg);
                         const __mmask8 m0 = _mm512_cmp_pd_mask(arg,zero,_CMP_LE_OQ);
                         const __mmask8 m1 = _mm512_cmp_pd_mask(arg,half,_CMP_LT_OQ);
                         const __mmask8 m2 = _mm512_cmp_pd_mask(
                                                          _mm512_mul_pd(ax,xinf),pi2,_CMP_LT_OQ);
                         const bool b      = m0 || (m1 && m2);
                         if(jint==1 && b) {
                             result = negate_zmm8r8(xinf);
                             return (result);
                         }                 
                         else if(_mm512_cmp_pd_mask(xmax,ax,_CMP_LT_OQ)) {
                             result = zero;
                             return (result);
                         }
                         if(_mm512_cmp_pd_mask(eight,ax,_CMP_LT_OQ)) {
                             goto L800;
                         }
                         else if(_mm512_cmp_pd_mask(ax,xsmall,_CMP_LE_OQ)) {
                                 if(jint==0) {
                                    result = _mm512_mul_pd(arg,half);
                                    return (result);
                                 }
                                 else {
                                    result = _mm512_div_pd(negate_zmm8r8(pi2),ax);
                                    return (result);
                                 }
                         }
                         
                         /*
                              Calculate J1 for appropriate interval, preserving
                              !  accuracy near the zero of J1.
                         */
                         
                         zsq = _mm512_mul_pd(ax,ax);
                         if(_mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ)) {
                             xnum = _mm512_fmadd_pd(pj0[6],zsq,
                                                       _mm512_fmadd_pd(pj0[5],zsq,pj0[4]));
                             xden = _mm512_add_pd(zsq,qj0[4]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[0]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[0]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[1]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[1]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[2]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[2]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[3]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[3]);
                             t0   = _mm512_sub_pd(ax,_mm512_div_pd(xj01,two56));
                             t1   = _mm512_add_pd(ax,xj0);
                             prod = _mm512_mul_pd(_mm512_sub_pd(t0,xj02),t1);
                         }
                         else {
                             xnum = pj1[0];
                             xden = _mm512_mul_pd(_mm512_add_pd(zsq,qj1[6],
                                                  _mm512_add_pd(zsq,qj1[0]));
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[1]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[1]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[2]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[2]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[3]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[3]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[4]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[4]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[5]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[5]);
                             t0   = _mm512_mul_pd(xnum,_mm512_sub_pd(ax,eight));
                             t1   = _mm512_add_pd(_mm512_add_pd(ax,eight),pj1[6]);
                             xnum = _mm512_mul_pd(t0,t1);
                             t0   = _mm512_mul_pd(xnum,_mm512_sub_pd(ax,four));
                             t1   = _mm512_add_pd(_mm512_add_pd(ax,four),pj1[7]);
                             xnum = _mm512_mul_pd(t0,t1);
                             t0   = _mm512_sub_pd(_mm512_sub_pd(ax,
                                                 _mm512_div_pd(xj11,two56)),xj12);
                             t1   = _mm512_add_pd(ax,xj1);
                             prod = _mm512_mul_pd(_mm512_mul_pd(arg,t0),t1);
                         }
                         result = _mm512_mul_pd(prod,_mm512_div_pd(xnum,xden));
                         if(jint==0) {
                             return (result);
                         }
                         
                         /*
                             Calculate Y1.  First find RESJ = pi/2 ln(x/xn) J1(x),
                             !  where xn is a zero of Y1.
                         */
                         
                          __mmask8 m = _mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ);
                          up = _mm512_mask_blend_pd(m,
                                               _mm512_sub_pd(_mm512_sub_pd(ax,
                                                                    _mm512_div_pd(xy01,two56)),xy02),
                                               _mm512_sub_pd(_mm512_sub_pd(ax,
                                                                    _mm512_div_pd(xy11,two56)),xy12));
                          xy = _mm512_mask_blend_pd(m,xy1,xy0);
                          
                          down = _mm512_add_pd(ax,xy);
                          if(_mm512_cmp_pd_mask(_mm512_abs_pd(up),
                                                _mm512_mul_pd(p17,down),
                                                                 _CMP_LT_OQ)) {
                              w   = _mm512_div_pd(up,down);
                              wsq = _mm512_mul_pd(w,w);
                              xnum= plg[0];
                              xden= _mm512_add_pd(wsq,qlg[0]);
                              xnum= _mm512_fmadd_pd(xnum,wsq,plg[1]);
                              xden= _mm512_fmadd_pd(xden,wsq,qlg[1]);
                              xnum= _mm512_fmadd_pd(xnum,wsq,plg[2]);
                              xden= _mm512_fmadd_pd(xden,wsq,qlg[2]);
                              xnum= _mm512_fmadd_pd(xnum,wsq,plg[3]);
                              xden= _mm512_fmקברטadd_pd(xden,wsq,qlg[3]); 
                              t0  = _mm512_mul_pd(w,_mm512_div_pd(xnum,xden));
                              t1  = _mm512_mul_pd(pi2,result);
                              resj= _mm512_mul_pd(t0,t1);                                    
                          }
                          else {
                              t0  = xlog(_mm512_div_pd(ax,xy));
                              resj= _mm512_mul_pd(pi2,
                                               _mm512_mul_pd(result,t0));
                          }
                          
                          /*
                             Now calculate Y1 for appropriate interval, preserving
                             !  accuracy near the zero of Y1.
                          */
                          
                          if(_mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ)) {
                              xnum = _mm512_fmadd_pd(py0[6],zsq,py0[0]);
                              xden = _mm512_add_pd(zsq,qy0[0]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[1]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[1]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[2]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[2]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[3]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[3]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[4]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[4]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[5]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[5]);
                          }
                          else {
                              xnum = _mm512_fmadd_pd(py1[8],zsq,py1[0]);
                              xden = _mm512_add_pd(zsq,qy1[0]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[1]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[1]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[2]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[2]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[3]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[3]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[4]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[4]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[5]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[5]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[6]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[6]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[7]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[7]);
                          }
                          t0       = _mm512_div_pd(xnum,xden);
                          t1       = _mm512_mul_pd(up,
                                                _mm512_div_pd(down,ax));
                          result   = _mm512_fmadd_pd(t0,t1,resj);
                          return (result);
                      L800:
                          z          = _mm5512_div_pd(eight,ax);
                          __m512i ti = mm512_cvt_roundpd_epi64(_mm512_div_pd(ax,twopi),
                                                     _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                          w          = _mm512_add_pd(_mm512_castsi512_pd(ti),throv8);
                          w          = _mm512_fmsub_pd(_mm512_sub_pd(ax,w),twopi1,
                                                                       _mm512_mul_pd(w,twopi2));
                          zsq        = _mm512_mul_pd(z,z);
                          xnum       = p0[5];
                          xden       = _mm512_add_pd(zsq,q0[5]);
                          up         = p1[5];
                          down       = _mm512_add_pd(zsq,q1[5]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[0]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[0]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[1]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[1]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[2]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[2]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[3]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[3]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[4]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[4]);
                          r0         = _mm512_div_pd(xnum,xden);
                          r1         = _mm512_div_pd(up,down);
                          t0         = _mm512_div_pd(rtpi2,_mm512_sqrt_pd(ax));
                          t1         = xsin(w);
                          register __m512d t2 = xcos(w);
                          register __m512d t3 = _mm512_mul_pd(z,r1);
                          if(jint==1) {
                             result = _mm512_mul_pd(t0,
                                               _mm12_fmsub_pd(r0,t2,
                                                            _mm512_mul_pd(t3,t1)));
                          }
                          else {
                             result = _mm512_mul_pd(t0,
                                               _mm12_fmsub_pd(r0,t1,
                                                            _mm512_mul_pd(t3,t2)));
                          }
                          if(jint==0 && 
                                _mm512_cmp+pd_mask(arg,zero,_CMP_LT_OQ)) {
                                 result = negate_zmm8r8(result);     
                          }
                          
                          return (result);
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d caljy1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) parg,
	                                 const int32_t jint) {
	                                 
	                 __ATTR_ALIGN__(64) const static 
	                 __m512d plg[4]  = {_mm512_set1_pd(-2.4562334077563243311e+1),
	                                    _mm512_set1_pd(2.3642701335621505212e+2),
                                            _mm512_set1_pd(-5.4989956895857911039e+2),
                                            _mm512_set1_pd(3.5687548468071500413e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qlg[4]  = {_mm512_set1_pd(-3.5553900764052419184e+1),
	                                    _mm512_set1_pd(1.9400230218539473193e+2),
                                            _mm512_set1_pd(-3.3442903192607538956e+2),
                                            _mm512_set1_pd(1.7843774234035750207e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj0[7]  = {_mm512_set1_pd(9.8062904098958257677e+5,
	                                    _mm512_set1_pd(-1.1548696764841276794e+8), 
                                            _mm512_set1_pd(6.6781041261492395835e+9),
                                            _mm512_set1_pd(-1.4258509801366645672e+11), 
                                            _mm512_set1_pd(-4.4615792982775076130e+3), 
                                            _mm512_set1_pd(1.0650724020080236441e+1),
                                            _mm512_set1_pd(-1.0767857011487300348e-2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj0[5]  = {_mm512_set1_pd(5.9117614494174794095e+5), 
	                                    _mm512_set1_pd(2.0228375140097033958e+8), 
                                            _mm512_set1_pd(4.2091902282580133541e+10), 
                                            _mm512_set1_pd(4.1868604460820175290e+12), 
                                            _mm512_set1_pd(1.0742272239517380498e+03)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d pj1[8]  = {_mm512_set1_pd(4.6179191852758252280e+00),
	                                    _mm512_set1_pd(-7.1329006872560947377e+3),
                                            _mm512_set1_pd(4.5039658105749078904e+6),
                                            _mm512_set1_pd(-1.4437717718363239107e+9),
                                            _mm512_set1_pd(2.3569285397217157313e+11),
                                            _mm512_set1_pd(-1.6324168293282543629e+13),
                                            _mm512_set1_pd(1.1357022719979468624e+14), 
                                            _mm512_set1_pd(1.0051899717115285432e+15)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qj1[7]  = {_mm512_set1_pd(1.1267125065029138050e+6), 
	                                    _mm512_set1_pd(6.4872502899596389593e+8),
                                            _mm512_set1_pd(2.7622777286244082666e+11), 
                                            _mm512_set1_pd(8.4899346165481429307e+13),
                                            _mm512_set1_pd(1.7128800897135812012e+16), 
                                            _mm512_set1_pd(1.7253905888447681194e+18), 
                                            _mm512_set1_pd(1.3886978985861357615e+3)}; 
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py0[7]  = {_mm512_set1_pd(2.2157953222280260820e+5),
	                                    _mm512_set1_pd(-5.9157479997408395984e+7), 
                                            _mm512_set1_pd(7.2144548214502560419e+9),
                                            _mm512_set1_pd(-3.7595974497819597599e+11),
                                            _mm512_set1_pd(5.4708611716525426053e+12), 
                                            _mm512_set1_pd(4.0535726612579544093e+13), 
                                            _mm512_set1_pd(-3.1714424660046133456e+2)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy0[6]  = {_mm512_set1_pd(8.2079908168393867438e+2), 
	                                    _mm512_set1_pd(3.8136470753052572164e+5),
                                            _mm512_set1_pd(1.2250435122182963220e+8), 
                                            _mm512_set1_pd(2.7800352738690585613e+10),
                                            _mm512_set1_pd(4.1272286200406461981e+12), 
                                            _mm512_set_pd(3.0737873921079286084e+14)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d py1[9]  =  {_mm512_set1_pd(1.9153806858264202986e+6),
	                                     _mm512_set1_pd(-1.1957961912070617006e+9), 
                                             _mm512_set1_pd(3.7453673962438488783e+11),
                                             _mm512_set1_pd(-5.9530713129741981618e+13), 
                                             _mm512_set1_pd(4.0686275289804744814e+15),
                                             _mm512_set1_pd(-2.3638408497043134724e+16),
                                             _mm512_set1_pd(-5.6808094574724204577e+18), 
                                             _mm512_set1_pd(1.1514276357909013326e+19), 
                                             _mm512_set1_pd(-1.2337180442012953128e+3)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d qy1[8]  =  {_mm512_set1_pd(1.2855164849321609336e+3), 
	                                     _mm512_set1_pd(1.0453748201934079734e+6), 
                                             _mm512_set1_pd(6.3550318087088919566e+8), 
                                             _mm512_set1_pd(3.0221766852960403645e+11), 
                                             _mm512_set1_pd(1.1187010065856971027e+14), 
                                             _mm512_set1_pd(3.0837179548112881950e+16),
                                             _mm512_set1_pd(5.6968198822857178911e+18), 
                                             _mm512_set1_pd(5.3321844313316185697e+20)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p0[6]  =   {_mm512_set1_pd(-1.0982405543459346727e+5),
	                                     _mm512_set1_pd(-1.5235293511811373833e+6),
                                             _mm512_set1_pd(-6.6033732483649391093e06),
                                             _mm512_set1_pd(-9.9422465050776411957e+6),
                                             _mm512_set1_pd(-4.4357578167941278571e+6),
                                             _mm512_set1_pd(-1.6116166443246101165e+3)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q0[6]  =   {_mm512_set1_pd(-1.0726385991103820119e+5),
	                                     _mm512_set1_pd(-1.5118095066341608816e+6),
                                             _mm512_set1_pd(-6.5853394797230870728e+6),
                                             _mm512_set1_pd(-9.9341243899345856590e+6), 
                                             _mm512_set1_pd(-4.4357578167941278568e+6),
                                             _mm512_set1_pd(-1.4550094401904961825e+3)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p1[6]  =   {_mm512_set1_pd(1.7063754290207680021e+3), 
	                                     _mm512_set1_pd(1.8494262873223866797e+4), 
                                             _mm512_set1_pd(6.6178836581270835179e+4), 
                                             _mm512_set1_pd(8.5145160675335701966e+4),
                                             _mm512_set1_pd(3.3220913409857223519e+4), 
                                             _mm512_set1_pd(3.5265133846636032186e+1)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q1[6]  =   {_mm512_set1_pd(3.7890229745772202641e+4), 
	                                     _mm512_set1_pd(4.0029443582266975117e+5),
                                             _mm512_set1_pd(1.4194606696037208929e+6), 
                                             _mm512_set1_pd(1.8194580422439972989e+6),
                                             _mm512_set1_pd(7.0871281941028743574e+5), 
                                             _mm512_set1_pd(8.6383677696049909675e+2)};
                         const __m512d eight = _mm512_set1_pd(8.0e+0);
                         const __m512d four  = _mm512_set1_pd(4.0e+0);
                         const __m512d  half = _mm512_set1_pd(0.5);
                         const __m512d throv8= _mm512_set1_pd(0.375);
                         const __m512d pi2   = _mm512_set1_pd(6.3661977236758134308e-1);
                         const __m512d p17   = _mm512_set1_pd(1.716e-1);
                         const __m512d twopi = _mm512_set1_pd(6.2831853071795864769e+0);
                         const __m512d zero  = _mm512_set1_pd(0.0);
                         const __m512d twopi1= _mm512_set1_pd(6.28125e+0);
                         const __m512d twopi2= _mm512_set1_pd(1.9353071795864769253e-3);
                         const __m512d two56 = _mm512_set1_pd(256.0e+0);
                         const __m512d rtpi2 = _mm512_set1_pd(7.9788456080286535588e-1);
                         const __m512d xmax  = _mm512_set1_pd(1.07e+9);
                         const __m512d xsmall= _mm512_set1_pd(9.31e-10);
                         const __m512d xinf  = _mm512_set1_pd(1.7e+38);  
                         const __m512d xj0   = _mm12_set1_pd(3.8317059702075123156e+0);
                         const __m512d xj1   = _mm512_set1_pd(7.0155866698156187535e+0);
                         const __m512d xy0   = _mm512_set1_pd(2.1971413260310170351e+0);
                         const __m512d xy1   = _mm512_set1_pd(5.4296810407941351328e+0);
                         const __m512d xj01  = _mm512_set1_pd(981.0e+0);
                         const __m512d xj02  = _mm512_set1_pd(-3.2527979248768438556e-4);
                         const __m512d xj11  = _mm512_set1_pd(1796.0e+0);
                         const __m512d xj12  = _mm512_set1_pd(-3.8330184381246462950e-5);
                         const __m512d xy01  = _mm512_set1_pd(562.0e+0);
                         const __m512d xy02  = _mm512_set1_pd(1.8288260310170351490e-3);
                         const __m512d xy11  = _mm512_set1_pd(1390.0e+0);
                         const __m512d xy12  = _mm512_set1_pd(-6.4592058648672279948e-6);
                         register __m512d arg,ax,down,prod,resj,result;
                         register __m512d r0,r1,up,w,wsq;
                         register __m512d xden,xnum,t0,t1,z,zsq;
                         
                         arg = _mm512_load_pd(&parg[0]);
                         ax = _mm512_abs_pd(arg);
                         const __mmask8 m0 = _mm512_cmp_pd_mask(arg,zero,_CMP_LE_OQ);
                         const __mmask8 m1 = _mm512_cmp_pd_mask(arg,half,_CMP_LT_OQ);
                         const __mmask8 m2 = _mm512_cmp_pd_mask(
                                                          _mm512_mul_pd(ax,xinf),pi2,_CMP_LT_OQ);
                         const bool b      = m0 || (m1 && m2);
                         if(jint==1 && b) {
                             result = negate_zmm8r8(xinf);
                             return (result);
                         }                 
                         else if(_mm512_cmp_pd_mask(xmax,ax,_CMP_LT_OQ)) {
                             result = zero;
                             return (result);
                         }
                         if(_mm512_cmp_pd_mask(eight,ax,_CMP_LT_OQ)) {
                             goto L800;
                         }
                         else if(_mm512_cmp_pd_mask(ax,xsmall,_CMP_LE_OQ)) {
                                 if(jint==0) {
                                    result = _mm512_mul_pd(arg,half);
                                    return (result);
                                 }
                                 else {
                                    result = _mm512_div_pd(negate_zmm8r8(pi2),ax);
                                    return (result);
                                 }
                         }
                         
                         /*
                              Calculate J1 for appropriate interval, preserving
                              !  accuracy near the zero of J1.
                         */
                         
                         zsq = _mm512_mul_pd(ax,ax);
                         if(_mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ)) {
                             xnum = _mm512_fmadd_pd(pj0[6],zsq,
                                                       _mm512_fmadd_pd(pj0[5],zsq,pj0[4]));
                             xden = _mm512_add_pd(zsq,qj0[4]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[0]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[0]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[1]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[1]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[2]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[2]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj0[3]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj0[3]);
                             t0   = _mm512_sub_pd(ax,_mm512_div_pd(xj01,two56));
                             t1   = _mm512_add_pd(ax,xj0);
                             prod = _mm512_mul_pd(_mm512_sub_pd(t0,xj02),t1);
                         }
                         else {
                             xnum = pj1[0];
                             xden = _mm512_mul_pd(_mm512_add_pd(zsq,qj1[6],
                                                  _mm512_add_pd(zsq,qj1[0]));
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[1]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[1]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[2]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[2]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[3]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[3]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[4]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[4]);
                             xnum = _mm512_fmadd_pd(xnum,zsq,pj1[5]);
                             xden = _mm512_fmadd_pd(xden,zsq,qj1[5]);
                             t0   = _mm512_mul_pd(xnum,_mm512_sub_pd(ax,eight));
                             t1   = _mm512_add_pd(_mm512_add_pd(ax,eight),pj1[6]);
                             xnum = _mm512_mul_pd(t0,t1);
                             t0   = _mm512_mul_pd(xnum,_mm512_sub_pd(ax,four));
                             t1   = _mm512_add_pd(_mm512_add_pd(ax,four),pj1[7]);
                             xnum = _mm512_mul_pd(t0,t1);
                             t0   = _mm512_sub_pd(_mm512_sub_pd(ax,
                                                 _mm512_div_pd(xj11,two56)),xj12);
                             t1   = _mm512_add_pd(ax,xj1);
                             prod = _mm512_mul_pd(_mm512_mul_pd(arg,t0),t1);
                         }
                         result = _mm512_mul_pd(prod,_mm512_div_pd(xnum,xden));
                         if(jint==0) {
                             return (result);
                         }
                         
                         /*
                             Calculate Y1.  First find RESJ = pi/2 ln(x/xn) J1(x),
                             !  where xn is a zero of Y1.
                         */
                         
                          __mmask8 m = _mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ);
                          up = _mm512_mask_blend_pd(m,
                                               _mm512_sub_pd(_mm512_sub_pd(ax,
                                                                    _mm512_div_pd(xy01,two56)),xy02),
                                               _mm512_sub_pd(_mm512_sub_pd(ax,
                                                                    _mm512_div_pd(xy11,two56)),xy12));
                          xy = _mm512_mask_blend_pd(m,xy1,xy0);
                          
                          down = _mm512_add_pd(ax,xy);
                          if(_mm512_cmp_pd_mask(_mm512_abs_pd(up),
                                                _mm512_mul_pd(p17,down),
                                                                 _CMP_LT_OQ)) {
                              w   = _mm512_div_pd(up,down);
                              wsq = _mm512_mul_pd(w,w);
                              xnum= plg[0];
                              xden= _mm512_add_pd(wsq,qlg[0]);
                              xnum= _mm512_fmadd_pd(xnum,wsq,plg[1]);
                              xden= _mm512_fmadd_pd(xden,wsq,qlg[1]);
                              xnum= _mm512_fmadd_pd(xnum,wsq,plg[2]);
                              xden= _mm512_fmadd_pd(xden,wsq,qlg[2]);
                              xnum= _mm512_fmadd_pd(xnum,wsq,plg[3]);
                              xden= _mm512_fmקברטadd_pd(xden,wsq,qlg[3]); 
                              t0  = _mm512_mul_pd(w,_mm512_div_pd(xnum,xden));
                              t1  = _mm512_mul_pd(pi2,result);
                              resj= _mm512_mul_pd(t0,t1);                                    
                          }
                          else {
                              t0  = xlog(_mm512_div_pd(ax,xy));
                              resj= _mm512_mul_pd(pi2,
                                               _mm512_mul_pd(result,t0));
                          }
                          
                          /*
                             Now calculate Y1 for appropriate interval, preserving
                             !  accuracy near the zero of Y1.
                          */
                          
                          if(_mm512_cmp_pd_mask(ax,four,_CMP_LE_OQ)) {
                              xnum = _mm512_fmadd_pd(py0[6],zsq,py0[0]);
                              xden = _mm512_add_pd(zsq,qy0[0]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[1]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[1]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[2]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[2]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[3]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[3]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[4]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[4]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py0[5]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy0[5]);
                          }
                          else {
                              xnum = _mm512_fmadd_pd(py1[8],zsq,py1[0]);
                              xden = _mm512_add_pd(zsq,qy1[0]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[1]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[1]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[2]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[2]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[3]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[3]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[4]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[4]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[5]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[5]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[6]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[6]);
                              xnum = _mm512_fmadd_pd(xnum,zsq,py1[7]);
                              xden = _mm512_fmadd_pd(xden,zsq,qy1[7]);
                          }
                          t0       = _mm512_div_pd(xnum,xden);
                          t1       = _mm512_mul_pd(up,
                                                _mm512_div_pd(down,ax));
                          result   = _mm512_fmadd_pd(t0,t1,resj);
                          return (result);
                      L800:
                          z          = _mm5512_div_pd(eight,ax);
                          __m512i ti = mm512_cvt_roundpd_epi64(_mm512_div_pd(ax,twopi),
                                                     _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
                          w          = _mm512_add_pd(_mm512_castsi512_pd(ti),throv8);
                          w          = _mm512_fmsub_pd(_mm512_sub_pd(ax,w),twopi1,
                                                                       _mm512_mul_pd(w,twopi2));
                          zsq        = _mm512_mul_pd(z,z);
                          xnum       = p0[5];
                          xden       = _mm512_add_pd(zsq,q0[5]);
                          up         = p1[5];
                          down       = _mm512_add_pd(zsq,q1[5]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[0]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[0]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[1]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[1]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[2]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[2]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[3]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[3]);
                          xnum       = _mm512_fmadd_pd(xnum,zsq,p0[4]);
                          xden       = _mm512_fmadd_pd(xden,zsq,q0[4]);
                          r0         = _mm512_div_pd(xnum,xden);
                          r1         = _mm512_div_pd(up,down);
                          t0         = _mm512_div_pd(rtpi2,_mm512_sqrt_pd(ax));
                          t1         = xsin(w);
                          register __m512d t2 = xcos(w);
                          register __m512d t3 = _mm512_mul_pd(z,r1);
                          if(jint==1) {
                             result = _mm512_mul_pd(t0,
                                               _mm12_fmsub_pd(r0,t2,
                                                            _mm512_mul_pd(t3,t1)));
                          }
                          else {
                             result = _mm512_mul_pd(t0,
                                               _mm12_fmsub_pd(r0,t1,
                                                            _mm512_mul_pd(t3,t2)));
                          }
                          if(jint==0 && 
                                _mm512_cmp+pd_mask(arg,zero,_CMP_LT_OQ)) {
                                 result = negate_zmm8r8(result);     
                          }
                          
                          return (result);
	         }
	         
	         
	         /*
	         !*****************************************************************************80
!
!! DAW evaluates Dawson's integral function.
!
!  Discussion:
!
!    This routine evaluates Dawson's integral,
!
!      F(x) = exp ( - x * x ) * Integral ( 0 <= t <= x ) exp ( t * t ) dt
!
!    for a real argument x.
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
!  Reference:
!
!    William Cody, Kathleen Paciorek, Henry Thacher,
!    Chebyshev Approximations for Dawson's Integral,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 171-178.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument of the function.
!
!    Output, real ( kind = 8 ) DAW, the value of the function.
	         */
	          
	            
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d dawson_zmm8r8(const __m512d xx) {
	           
	                 __ATTR_ALIGN__(64) const static 
	                 __m512d p1[10]  = {_mm512_set1_pd(-2.69020398788704782410e-12), 
	                                    _mm512_set1_pd(4.18572065374337710778e-10),
                                            _mm512_set1_pd(-1.34848304455939419963e-8), 
                                            _mm512_set1_pd(9.28264872583444852976e-7),
                                            _mm512_set1_pd(-1.23877783329049120592e-5), 
                                            _mm512_set1_pd(4.07205792429155826266e-4),
                                            _mm512_set1_pd(-2.84388121441008500446e-3), 
                                            _mm512_set1_pd(4.70139022887204722217e-2), 
                                            _mm512_set1_pd(-1.38868086253931995101e-1), 
                                            _mm512_set1_pd(1.00000000000000000004e+00)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q1[10]  = {_mm512_set1_pd(1.71257170854690554214e-10), 
	                                    _mm512_set1_pd(1.19266846372297253797e-8), 
                                            _mm512_set1_pd(4.32287827678631772231e-7), 
                                            _mm12_set1_pd(1.03867633767414421898e-5), 
                                            _mm5512_set1_pd(1.78910965284246249340e-4), 
                                            _mm512_set1_pd(2.26061077235076703171e-3), 
                                            _mm512_set1_pd(2.07422774641447644725e-2), 
                                            _mm512_set1_pd(1.32212955897210128811e-1), 
                                            _mm512_set1_pd(5.27798580412734677256e-1), 
                                            _mm512_set1_pd(1.00000000000000000000e+00)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p2[10]  = {_mm512_set1_pd(-1.70953804700855494930e+0),
	                                    _mm512_set1_pd(-3.79258977271042880786e+1), 
                                            _mm512_set1_pd(2.61935631268825992835e+1), 
                                            _mm512_set1_pd(1.25808703738951251885e+1), 
                                            _mm512_set1_pd(-2.27571829525075891337e+1), 
                                            _mm512_set1_pd(4.56604250725163310122e+00), 
                                            _mm512_set1_pd(-7.33080089896402870750e+00), 
                                            _mm512_set1_pd(4.65842087940015295573e+01),
                                            _mm512_set1_pd(-1.73717177843672791149e+01), 
                                            _mm512_set1_pd(5.00260183622027967838e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q2[9]  = {_mm512_set1_pd(1.82180093313514478378e+00), 
	                                   _mm512_set1_pd(1.10067081034515532891e+03), 
                                           _mm512_set1_pd(-7.08465686676573000364e+00), 
                                           _mm512_set1_pd(4.53642111102577727153e+02), 
                                           _mm512_set1_pd(4.06209742218935689922e+01), 
                                           _mm512_set1_pd(3.02890110610122663923e+02),
                                           _mm512_set1_pd(1.70641269745236227356e+02), 
                                           _mm512_set1_pd(9.51190923960381458747e+02),
                                           _mm512_set1_pd(2.06522691539642105009e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p3[10]  = {_mm512_set1_pd(4.55169503255094815112e+00),
	                                    _mm512_set1_pd(-1.86647123338493852582e+01),
                                            _mm512_set1_pd(-7.36315669126830526754e+00),
                                            _mm512_set1_pd(-6.68407240337696756838e+01), 
                                            _mm512_set1_pd(4.84507265081491452130e+01), 
                                            _mm512_set1_pd(2.69790586735467649969e+01),
                                            _mm512_set1_pd(-3.35044149820592449072e+01), 
                                            _mm512_set1_pd(7.50964459838919612289e+00),
                                            _mm512_set1_pd(-1.48432341823343965307e+00), 
                                            _mm512_set1_pd(4.99999810924858824981e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q3[9]  = {_mm512_set1_pd(4.47820908025971749852e+01), 
	                                   _mm512_set1_pd(9.98607198039452081913e+01),
                                           _mm512_set1_pd(1.40238373126149385228e+01), 
                                           _mm512_set1_pd(3.48817758822286353588e+03),
                                           _mm512_set1_pd(-9.18871385293215873406e+00), 
                                           _mm512_set1_pd(1.24018500009917163023e+03),
                                           _mm512_set1_pd(-6.88024952504512254535e+01),
                                           _mm512_set1_pd(-2.31251575385145143070e+00),
                                           _mm512_set1_pd(2.50041492369922381761e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p4[10]  = {_mm512_set1_pd(-8.11753647558432685797e+00),
	                                    _mm512_set1_pd(-3.84043882477454453430e+01),
                                            _mm512_set1_pd(-2.23787669028751886675e+01),
                                            _mm512_set1_pd(-2.88301992467056105854e+01),
                                            _mm512_set1_pd(-5.99085540418222002197e+00),
                                            _mm512_set1_pd(-1.13867365736066102577e+01), 
                                            _mm512_set1_pd(-6.52828727526980741590e+00),
                                            _mm512_set1_pd(-4.50002293000355585708e+00),
                                            _mm512_set1_pd(-2.50000000088955834952e+00), 
                                            _mm512_set1_pd(5.00000000000000488400e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q4[9]  = {_mm512_set1_pd(2.69382300417238816428e+02), 
	                                   _mm512_set1_pd(5.04198958742465752861e+01),
                                           _mm512_set1_pd(6.11539671480115846173e+01), 
                                           _mm512_set1_pd(2.08210246935564547889e+02),
                                           _mm512_set1_pd(1.97325365692316183531e+01),
                                           _mm512_set1_pd(-1.22097010558934838708e+01),
                                           _mm512_set1_pd(-6.99732735041547247161e+00),
                                           _mm512_set1_pd(-2.49999970104184464568e+00), 
                                           _mm512_set1_pd(7.49999999999027092188e-01)};
                         const __m512d zero  = _mm512_set1_pd(0.0e+00);
                         const __m512d half  = _mm5512_set1_pd(0.5e+00);
                         const __m512d one   = _mm512_set1_pd(1.0e+00);
                         const __m512d six25 = _mm512_set1_pd(6.25e+00);
                         const __m512d one225= _mm512_set1_pd(12.25e+0);
                         const __m512d two5  = _mm512_set1_pd(25.0e+0);
                         const __m512d xsmall= _mm512_set1_pd(1.05e-08);
                         const __m512d xlarge= _mm512_set1_pd(9.49e+07);
                         const __m512d xmax  = _mm512_set1_pd(2.24e+307);
                         register __m512d daw,frac,sump,sumq,w2,x,y,t0,t1;
                        
                         x = xx;
                         t0= _mm512_abs_pd(x);
                         if(_mm512_cmp_pd_mask(xlarge,t0,_CMP_LT_OQ)) {
                             const __mmask8 m = _mm512_cmp_pd_mask(t0,xmax,_CMP_LE_OQ);
                             daw = _mm512_mask_blend_pd(m,zero,
                                                   _mm512_div_pd(half,x));         
                         } 
                         else if(_mm512_cmp_pd_mask(t0,xsmall,_CMP_LT_OQ)) {
                             daw = x;
                         }
                         else {
                             y = _mm512_mul_pd(x,x);
                             if(_mm512_cmp_pd_mask(y,six25,_CMP_LT_OQ)) {
                                sump = p1[0];
                                sumq = q1[0];
                                sump = _mm512_fmadd_pd(sump,y,p1[1]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[1]);
                                sump = _mm512_fmadd_pd(sump,y,p1[2]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[2]);
                                sump = _mm512_fmadd_pd(sump,y,p1[3]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[3]);
                                sump = _mm512_fmadd_pd(sump,y,p1[4]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[4]);
                                sump = _mm512_fmadd_pd(sump,y,p1[5]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[5]);
                                sump = _mm512_fmadd_pd(sump,y,p1[6]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[6]);
                                sump = _mm512_fmadd_pd(sump,y,p1[7]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[7]);
                                sump = _mm512_fmadd_pd(sump,y,p1[8]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[8]);
                                sump = _mm512_fmadd_pd(sump,y,p1[9]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[9]);
                                daw  = _mm512_mul_pd(x,_mm512_div_pd(sump,sumq));
                             }
                             else if(_mm512_cmp_pd_mask(y,one225,_CMP_LT_OQ)) {
                                frac = zero;
                                t0   = _mm512_add_pd(y,frac);
                                frac = _mm512_div_pd(q2[0],
                                                 _mm512_add_pd(p2[0],t0));
                                frac = _mm512_div_pd(q2[1],
                                                 _mm512_add_pd(p2[1],t0));
                                frac = _mm512_div_pd(q2[2],
                                                 _mm512_add_pd(p2[2],t0));
                                frac = _mm512_div_pd(q2[3],
                                                 _mm512_add_pd(p2[3],t0));
                                frac = _mm512_div_pd(q2[4],
                                                 _mm512_add_pd(p2[4],t0));
                                frac = _mm512_div_pd(q2[5],
                                                 _mm512_add_pd(p2[5],t0));
                                frac = _mm512_div_pd(q2[6],
                                                 _mm512_add_pd(p2[6],t0));
                                frac = _mm512_div_pd(q2[7],
                                                 _mm512_add_pd(p2[7],t0));
                                frac = _mm512_div_pd(q2[8],
                                                 _mm512_add_pd(p2[8],t0));
                                daw  = _mm512_div_pd(_mm512_add_pd(p2[9],frac),x);
                             }   
                             else if(_mm512_cmp_pd_mask(y,two5,_CMP_LT_OQ)) {
                                frac = zero;
                                t0   = _mm512_add_pd(y,frac);
                                frac = _mm512_div_pd(q3[0],
                                                 _mm512_add_pd(p3[0],t0));
                                frac = _mm512_div_pd(q3[1],
                                                 _mm512_add_pd(p3[1],t0));
                                frac = _mm512_div_pd(q3[2],
                                                 _mm512_add_pd(p3[2],t0));
                                frac = _mm512_div_pd(q3[3],
                                                 _mm512_add_pd(p3[3],t0));
                                frac = _mm512_div_pd(q3[4],
                                                 _mm512_add_pd(p3[4],t0));
                                frac = _mm512_div_pd(q3[5],
                                                 _mm512_add_pd(p3[5],t0));
                                frac = _mm512_div_pd(q3[6],
                                                 _mm512_add_pd(p3[6],t0));
                                frac = _mm512_div_pd(q3[7],
                                                 _mm512_add_pd(p3[7],t0));
                                frac = _mm512_div_pd(q3[8],
                                                 _mm512_add_pd(p3[8],t0));
                                daw  = _mm512_div_pd(_mm512_add_pd(p3[9],frac),x);
                             } 
                             else {
                                w2   = _mm512_div_pd(_mm512_div_pd(one,x),x);
                                frac = zero;
                                t0   = _mm512_add_pd(y,frac);
                                frac = _mm512_div_pd(q4[0],
                                                 _mm512_add_pd(p4[0],t0));
                                frac = _mm512_div_pd(q4[1],
                                                 _mm512_add_pd(p4[1],t0));
                                frac = _mm512_div_pd(q4[2],
                                                 _mm512_add_pd(p4[2],t0));
                                frac = _mm512_div_pd(q4[3],
                                                 _mm512_add_pd(p4[3],t0));
                                frac = _mm512_div_pd(q4[4],
                                                 _mm512_add_pd(p4[4],t0));
                                frac = _mm512_div_pd(q4[5],
                                                 _mm512_add_pd(p4[5],t0));
                                frac = _mm512_div_pd(q4[6],
                                                 _mm512_add_pd(p4[6],t0));
                                frac = _mm512_div_pd(q4[7],
                                                 _mm512_add_pd(p4[7],t0));
                                frac = _mm512_div_pd(q4[8],
                                                 _mm512_add_pd(p4[8],t0));
                                frac = _mm512_add_pd(p4[9],frac);
                                t0   = _mm512_mul_pd(w2,frac);
                                daw  = _mm512_div_pd(_mm512_fmadd_pd(t0,half,half),x);
                             }                        
                         }
                         
                         return (daw);
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d dawson_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pxx) {
	           
	                 __ATTR_ALIGN__(64) const static 
	                 __m512d p1[10]  = {_mm512_set1_pd(-2.69020398788704782410e-12), 
	                                    _mm512_set1_pd(4.18572065374337710778e-10),
                                            _mm512_set1_pd(-1.34848304455939419963e-8), 
                                            _mm512_set1_pd(9.28264872583444852976e-7),
                                            _mm512_set1_pd(-1.23877783329049120592e-5), 
                                            _mm512_set1_pd(4.07205792429155826266e-4),
                                            _mm512_set1_pd(-2.84388121441008500446e-3), 
                                            _mm512_set1_pd(4.70139022887204722217e-2), 
                                            _mm512_set1_pd(-1.38868086253931995101e-1), 
                                            _mm512_set1_pd(1.00000000000000000004e+00)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q1[10]  = {_mm512_set1_pd(1.71257170854690554214e-10), 
	                                    _mm512_set1_pd(1.19266846372297253797e-8), 
                                            _mm512_set1_pd(4.32287827678631772231e-7), 
                                            _mm12_set1_pd(1.03867633767414421898e-5), 
                                            _mm5512_set1_pd(1.78910965284246249340e-4), 
                                            _mm512_set1_pd(2.26061077235076703171e-3), 
                                            _mm512_set1_pd(2.07422774641447644725e-2), 
                                            _mm512_set1_pd(1.32212955897210128811e-1), 
                                            _mm512_set1_pd(5.27798580412734677256e-1), 
                                            _mm512_set1_pd(1.00000000000000000000e+00)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p2[10]  = {_mm512_set1_pd(-1.70953804700855494930e+0),
	                                    _mm512_set1_pd(-3.79258977271042880786e+1), 
                                            _mm512_set1_pd(2.61935631268825992835e+1), 
                                            _mm512_set1_pd(1.25808703738951251885e+1), 
                                            _mm512_set1_pd(-2.27571829525075891337e+1), 
                                            _mm512_set1_pd(4.56604250725163310122e+00), 
                                            _mm512_set1_pd(-7.33080089896402870750e+00), 
                                            _mm512_set1_pd(4.65842087940015295573e+01),
                                            _mm512_set1_pd(-1.73717177843672791149e+01), 
                                            _mm512_set1_pd(5.00260183622027967838e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q2[9]  = {_mm512_set1_pd(1.82180093313514478378e+00), 
	                                   _mm512_set1_pd(1.10067081034515532891e+03), 
                                           _mm512_set1_pd(-7.08465686676573000364e+00), 
                                           _mm512_set1_pd(4.53642111102577727153e+02), 
                                           _mm512_set1_pd(4.06209742218935689922e+01), 
                                           _mm512_set1_pd(3.02890110610122663923e+02),
                                           _mm512_set1_pd(1.70641269745236227356e+02), 
                                           _mm512_set1_pd(9.51190923960381458747e+02),
                                           _mm512_set1_pd(2.06522691539642105009e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p3[10]  = {_mm512_set1_pd(4.55169503255094815112e+00),
	                                    _mm512_set1_pd(-1.86647123338493852582e+01),
                                            _mm512_set1_pd(-7.36315669126830526754e+00),
                                            _mm512_set1_pd(-6.68407240337696756838e+01), 
                                            _mm512_set1_pd(4.84507265081491452130e+01), 
                                            _mm512_set1_pd(2.69790586735467649969e+01),
                                            _mm512_set1_pd(-3.35044149820592449072e+01), 
                                            _mm512_set1_pd(7.50964459838919612289e+00),
                                            _mm512_set1_pd(-1.48432341823343965307e+00), 
                                            _mm512_set1_pd(4.99999810924858824981e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q3[9]  = {_mm512_set1_pd(4.47820908025971749852e+01), 
	                                   _mm512_set1_pd(9.98607198039452081913e+01),
                                           _mm512_set1_pd(1.40238373126149385228e+01), 
                                           _mm512_set1_pd(3.48817758822286353588e+03),
                                           _mm512_set1_pd(-9.18871385293215873406e+00), 
                                           _mm512_set1_pd(1.24018500009917163023e+03),
                                           _mm512_set1_pd(-6.88024952504512254535e+01),
                                           _mm512_set1_pd(-2.31251575385145143070e+00),
                                           _mm512_set1_pd(2.50041492369922381761e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d p4[10]  = {_mm512_set1_pd(-8.11753647558432685797e+00),
	                                    _mm512_set1_pd(-3.84043882477454453430e+01),
                                            _mm512_set1_pd(-2.23787669028751886675e+01),
                                            _mm512_set1_pd(-2.88301992467056105854e+01),
                                            _mm512_set1_pd(-5.99085540418222002197e+00),
                                            _mm512_set1_pd(-1.13867365736066102577e+01), 
                                            _mm512_set1_pd(-6.52828727526980741590e+00),
                                            _mm512_set1_pd(-4.50002293000355585708e+00),
                                            _mm512_set1_pd(-2.50000000088955834952e+00), 
                                            _mm512_set1_pd(5.00000000000000488400e-01)};
                         __ATTR_ALIGN__(64) const static 
	                 __m512d q4[9]  = {_mm512_set1_pd(2.69382300417238816428e+02), 
	                                   _mm512_set1_pd(5.04198958742465752861e+01),
                                           _mm512_set1_pd(6.11539671480115846173e+01), 
                                           _mm512_set1_pd(2.08210246935564547889e+02),
                                           _mm512_set1_pd(1.97325365692316183531e+01),
                                           _mm512_set1_pd(-1.22097010558934838708e+01),
                                           _mm512_set1_pd(-6.99732735041547247161e+00),
                                           _mm512_set1_pd(-2.49999970104184464568e+00), 
                                           _mm512_set1_pd(7.49999999999027092188e-01)};
                         const __m512d zero  = _mm512_set1_pd(0.0e+00);
                         const __m512d half  = _mm5512_set1_pd(0.5e+00);
                         const __m512d one   = _mm512_set1_pd(1.0e+00);
                         const __m512d six25 = _mm512_set1_pd(6.25e+00);
                         const __m512d one225= _mm512_set1_pd(12.25e+0);
                         const __m512d two5  = _mm512_set1_pd(25.0e+0);
                         const __m512d xsmall= _mm512_set1_pd(1.05e-08);
                         const __m512d xlarge= _mm512_set1_pd(9.49e+07);
                         const __m512d xmax  = _mm512_set1_pd(2.24e+307);
                         register __m512d xx,daw,frac,sump,sumq,w2,x,y,t0,t1;
                        
                         xx = _mm512_load_pd(&pxx[0]);
                         x = xx;
                         t0= _mm512_abs_pd(x);
                         if(_mm512_cmp_pd_mask(xlarge,t0,_CMP_LT_OQ)) {
                             const __mmask8 m = _mm512_cmp_pd_mask(t0,xmax,_CMP_LE_OQ);
                             daw = _mm512_mask_blend_pd(m,zero,
                                                   _mm512_div_pd(half,x));         
                         } 
                         else if(_mm512_cmp_pd_mask(t0,xsmall,_CMP_LT_OQ)) {
                             daw = x;
                         }
                         else {
                             y = _mm512_mul_pd(x,x);
                             if(_mm512_cmp_pd_mask(y,six25,_CMP_LT_OQ)) {
                                sump = p1[0];
                                sumq = q1[0];
                                sump = _mm512_fmadd_pd(sump,y,p1[1]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[1]);
                                sump = _mm512_fmadd_pd(sump,y,p1[2]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[2]);
                                sump = _mm512_fmadd_pd(sump,y,p1[3]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[3]);
                                sump = _mm512_fmadd_pd(sump,y,p1[4]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[4]);
                                sump = _mm512_fmadd_pd(sump,y,p1[5]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[5]);
                                sump = _mm512_fmadd_pd(sump,y,p1[6]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[6]);
                                sump = _mm512_fmadd_pd(sump,y,p1[7]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[7]);
                                sump = _mm512_fmadd_pd(sump,y,p1[8]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[8]);
                                sump = _mm512_fmadd_pd(sump,y,p1[9]);
                                sumq = _mm512_fmadd_pd(sumq,y,q1[9]);
                                daw  = _mm512_mul_pd(x,_mm512_div_pd(sump,sumq));
                             }
                             else if(_mm512_cmp_pd_mask(y,one225,_CMP_LT_OQ)) {
                                frac = zero;
                                t0   = _mm512_add_pd(y,frac);
                                frac = _mm512_div_pd(q2[0],
                                                 _mm512_add_pd(p2[0],t0));
                                frac = _mm512_div_pd(q2[1],
                                                 _mm512_add_pd(p2[1],t0));
                                frac = _mm512_div_pd(q2[2],
                                                 _mm512_add_pd(p2[2],t0));
                                frac = _mm512_div_pd(q2[3],
                                                 _mm512_add_pd(p2[3],t0));
                                frac = _mm512_div_pd(q2[4],
                                                 _mm512_add_pd(p2[4],t0));
                                frac = _mm512_div_pd(q2[5],
                                                 _mm512_add_pd(p2[5],t0));
                                frac = _mm512_div_pd(q2[6],
                                                 _mm512_add_pd(p2[6],t0));
                                frac = _mm512_div_pd(q2[7],
                                                 _mm512_add_pd(p2[7],t0));
                                frac = _mm512_div_pd(q2[8],
                                                 _mm512_add_pd(p2[8],t0));
                                daw  = _mm512_div_pd(_mm512_add_pd(p2[9],frac),x);
                             }   
                             else if(_mm512_cmp_pd_mask(y,two5,_CMP_LT_OQ)) {
                                frac = zero;
                                t0   = _mm512_add_pd(y,frac);
                                frac = _mm512_div_pd(q3[0],
                                                 _mm512_add_pd(p3[0],t0));
                                frac = _mm512_div_pd(q3[1],
                                                 _mm512_add_pd(p3[1],t0));
                                frac = _mm512_div_pd(q3[2],
                                                 _mm512_add_pd(p3[2],t0));
                                frac = _mm512_div_pd(q3[3],
                                                 _mm512_add_pd(p3[3],t0));
                                frac = _mm512_div_pd(q3[4],
                                                 _mm512_add_pd(p3[4],t0));
                                frac = _mm512_div_pd(q3[5],
                                                 _mm512_add_pd(p3[5],t0));
                                frac = _mm512_div_pd(q3[6],
                                                 _mm512_add_pd(p3[6],t0));
                                frac = _mm512_div_pd(q3[7],
                                                 _mm512_add_pd(p3[7],t0));
                                frac = _mm512_div_pd(q3[8],
                                                 _mm512_add_pd(p3[8],t0));
                                daw  = _mm512_div_pd(_mm512_add_pd(p3[9],frac),x);
                             } 
                             else {
                                w2   = _mm512_div_pd(_mm512_div_pd(one,x),x);
                                frac = zero;
                                t0   = _mm512_add_pd(y,frac);
                                frac = _mm512_div_pd(q4[0],
                                                 _mm512_add_pd(p4[0],t0));
                                frac = _mm512_div_pd(q4[1],
                                                 _mm512_add_pd(p4[1],t0));
                                frac = _mm512_div_pd(q4[2],
                                                 _mm512_add_pd(p4[2],t0));
                                frac = _mm512_div_pd(q4[3],
                                                 _mm512_add_pd(p4[3],t0));
                                frac = _mm512_div_pd(q4[4],
                                                 _mm512_add_pd(p4[4],t0));
                                frac = _mm512_div_pd(q4[5],
                                                 _mm512_add_pd(p4[5],t0));
                                frac = _mm512_div_pd(q4[6],
                                                 _mm512_add_pd(p4[6],t0));
                                frac = _mm512_div_pd(q4[7],
                                                 _mm512_add_pd(p4[7],t0));
                                frac = _mm512_div_pd(q4[8],
                                                 _mm512_add_pd(p4[8],t0));
                                frac = _mm512_add_pd(p4[9],frac);
                                t0   = _mm512_mul_pd(w2,frac);
                                daw  = _mm512_div_pd(_mm512_fmadd_pd(t0,half,half),x);
                             }                        
                         }
                         
                         return (daw);
	         }
	         
	        
	        /*
	        !*****************************************************************************80
!
!! CERF computes the error function and derivative for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ), the argument.
!
!    Output, complex ( kind = 8 ) CER, CDER, the values of erf(z) and erf'(z).
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cerf_zmm8r8(const __m512d zr,
	                            const __m512d zi,
	                            __m512d * __restrict cerr,
	                            __m512d * __restrict ceri,
	                            __m512d * __restrict cderr,
	                            __m512d * __restrict cderi) {
	                     
	                __ATTR_ALIGN__(64) const static
	                __m512d LUT1[100] = {_mm512_set1_pd(1.50000000),
	                                    _mm512_set1_pd(2.50000000),
	                                    _mm512_set1_pd(3.50000000),
	                                    _mm512_set1_pd(4.50000000),
	                                    _mm512_set1_pd(5.50000000),
	                                    _mm512_set1_pd(6.50000000),
	                                    _mm512_set1_pd(7.50000000),    
                                            _mm512_set1_pd(8.50000000),    
                                            _mm512_set1_pd(9.50000000),    
                                            _mm512_set1_pd(10.5000000),    
                                            _mm512_set1_pd(11.5000000),    
                                            _mm512_set1_pd(12.5000000),    
                                            _mm512_set1_pd(13.5000000),    
                                            _mm512_set1_pd(14.5000000),   
                                            _mm512_set1_pd(15.5000000),    
                                            _mm512_set1_pd(16.5000000),    
                                            _mm512_set1_pd(17.5000000),    
                                            _mm512_set1_pd(18.5000000),    
                                            _mm512_set1_pd(19.5000000),    
                                            _mm512_set1_pd(20.5000000),    
                                            _mm512_set1_pd(21.5000000),    
                                            _mm512_set1_pd(22.5000000),    
                                            _mm512_set1_pd(23.5000000),   
                                            _mm512_set1_pd(24.5000000),    
                                            _mm512_set1_pd(25.5000000),   
                                            _mm512_set1_pd(26.5000000),    
                                            _mm512_set1_pd(27.5000000),    
                                            _mm512_set1_pd(28.5000000),    
                                            _mm512_set1_pd(29.5000000),    
                                            _mm512_set1_pd(30.5000000),   
                                            _mm512_set1_pd(31.5000000),    
                                            _mm512_set1_pd(32.5000000),    
                                            _mm512_set1_pd(33.5000000),    
                                            _mm512_set1_pd(34.5000000),    
                                            _mm512_set1_pd(35.5000000),   
                                            _mm512_set1_pd(36.5000000),    
                                            _mm512_set1_pd(37.5000000),    
                                            _mm512_set1_pd(38.5000000),    
                                            _mm512_set1_pd(39.5000000),
                                            _mm512_set1_pd(40.5000000),    
                                            _mm512_set1_pd(41.5000000),    
                                            _mm512_set1_pd(42.5000000),    
                                            _mm512_set1_pd(43.5000000),    
                                            _mm512_set1_pd(44.5000000),    
                                            _mm512_set1_pd(45.5000000),    
                                            _mm512_set1_pd(46.5000000),    
                                            _mm512_set1_pd(47.5000000),    
                                            _mm512_set1_pd(48.5000000),    
                                            _mm512_set1_pd(49.5000000),    
                                            _mm512_set1_pd(50.5000000),    
                                            _mm512_set1_pd(51.5000000),    
                                            _mm512_set1_pd(52.5000000),    
                                            _mm512_set1_pd(53.5000000),    
                                            _mm512_set1_pd(54.5000000),    
                                            _mm512_set1_pd(55.5000000),    
                                            _mm512_set1_pd(56.5000000),   
                                            _mm512_set1_pd(57.5000000),    
                                            _mm512_set1_pd(58.5000000),    
                                            _mm512_set1_pd(59.5000000),    
                                            _mm512_set1_pd(60.5000000),    
                                            _mm512_set1_pd(61.5000000),    
                                            _mm512_set1_pd(62.5000000),    
                                            _mm512_set1_pd(63.5000000),    
                                            _mm512_set1_pd(64.5000000),    
                                            _mm512_set1_pd(65.5000000),    
                                            _mm512_set1_pd(66.5000000),    
                                            _mm512_set1_pd(67.5000000),    
                                            _mm512_set1_pd(68.5000000),    
                                            _mm512_set1_pd(69.5000000),    
                                            _mm512_set1_pd(70.5000000),    
                                            _mm512_set1_pd(71.5000000),    
                                            _mm512_set1_pd(72.5000000),    
                                            _mm512_set1_pd(73.5000000),    
                                            _mm512_set1_pd(74.5000000),    
                                            _mm512_set1_pd(75.5000000),    
                                            _mm512_set1_pd(76.5000000),    
                                            _mm512_set1_pd(77.5000000),    
                                            _mm512_set1_pd(78.5000000),    
                                            _mm512_set1_pd(79.5000000),    
                                            _mm512_set1_pd(80.5000000),    
                                            _mm512_set1_pd(81.5000000),    
                                            _mm512_set1_pd(82.5000000),    
                                            _mm512_set1_pd(83.5000000),    
                                            _mm512_set1_pd(84.5000000),    
                                            _mm512_set1_pd(85.5000000),    
                                            _mm512_set1_pd(86.5000000),    
                                            _mm512_set1_pd(87.5000000),    
                                            _mm512_set1_pd(88.5000000),    
                                            _mm512_set1_pd(89.5000000),    
                                            _mm512_set1_pd(90.5000000),    
                                            _mm512_set1_pd(91.5000000),    
                                            _mm512_set1_pd(92.5000000),    
                                            _mm512_set1_pd(93.5000000),    
                                            _mm512_set1_pd(94.5000000),    
                                            _mm512_set1_pd(95.5000000),    
                                            _mm512_set1_pd(96.5000000),    
                                            _mm512_set1_pd(97.5000000),    
                                            _mm512_set1_pd(98.5000000),    
                                            _mm512_set1_pd(99.5000000),    
                                            _mm512_set1_pd(100.500000)}; 
                        __ATTR_ALIGN__(64) const static
	                __m512d LUT2[100] = {_mm512_set1_pd(1.0000000),
	                                    _mm512_set1_pd(2.0000000),
	                                    _mm512_set1_pd(3.0000000),
	                                    _mm512_set1_pd(4.0000000),
	                                    _mm512_set1_pd(5.0000000),
	                                    _mm512_set1_pd(6.0000000),
	                                    _mm512_set1_pd(7.0000000),    
                                            _mm512_set1_pd(8.0000000),    
                                            _mm512_set1_pd(9.0000000),    
                                            _mm512_set1_pd(10.000000),    
                                            _mm512_set1_pd(11.000000),    
                                            _mm512_set1_pd(12.000000),    
                                            _mm512_set1_pd(13.000000),    
                                            _mm512_set1_pd(14.000000),   
                                            _mm512_set1_pd(15.000000),    
                                            _mm512_set1_pd(16.000000),    
                                            _mm512_set1_pd(17.000000),    
                                            _mm512_set1_pd(18.000000),    
                                            _mm512_set1_pd(19.000000),    
                                            _mm512_set1_pd(20.000000),    
                                            _mm512_set1_pd(21.000000),    
                                            _mm512_set1_pd(22.000000),    
                                            _mm512_set1_pd(23.000000),   
                                            _mm512_set1_pd(24.000000),    
                                            _mm512_set1_pd(25.000000),   
                                            _mm512_set1_pd(26.000000),    
                                            _mm512_set1_pd(27.000000),    
                                            _mm512_set1_pd(28.000000),    
                                            _mm512_set1_pd(29.000000),    
                                            _mm512_set1_pd(30.000000),   
                                            _mm512_set1_pd(31.000000),    
                                            _mm512_set1_pd(32.000000),    
                                            _mm512_set1_pd(33.000000),    
                                            _mm512_set1_pd(34.000000),    
                                            _mm512_set1_pd(35.000000),   
                                            _mm512_set1_pd(36.000000),    
                                            _mm512_set1_pd(37.000000),    
                                            _mm512_set1_pd(38.000000),    
                                            _mm512_set1_pd(39.000000),
                                            _mm512_set1_pd(40.000000),    
                                            _mm512_set1_pd(41.000000),    
                                            _mm512_set1_pd(42.000000),    
                                            _mm512_set1_pd(43.000000),    
                                            _mm512_set1_pd(44.000000),    
                                            _mm512_set1_pd(45.000000),    
                                            _mm512_set1_pd(46.000000),    
                                            _mm512_set1_pd(47.000000),    
                                            _mm512_set1_pd(48.000000),    
                                            _mm512_set1_pd(49.000000),    
                                            _mm512_set1_pd(50.000000),    
                                            _mm512_set1_pd(51.000000),    
                                            _mm512_set1_pd(52.000000),    
                                            _mm512_set1_pd(53.000000),    
                                            _mm512_set1_pd(54.000000),    
                                            _mm512_set1_pd(55.000000),    
                                            _mm512_set1_pd(56.000000),   
                                            _mm512_set1_pd(57.000000),    
                                            _mm512_set1_pd(58.000000),    
                                            _mm512_set1_pd(59.000000),    
                                            _mm512_set1_pd(60.000000),    
                                            _mm512_set1_pd(61.000000),    
                                            _mm512_set1_pd(62.000000),    
                                            _mm512_set1_pd(63.000000),    
                                            _mm512_set1_pd(64.000000),    
                                            _mm512_set1_pd(65.000000),    
                                            _mm512_set1_pd(66.000000),    
                                            _mm512_set1_pd(67.000000),    
                                            _mm512_set1_pd(68.000000),    
                                            _mm512_set1_pd(69.000000),    
                                            _mm512_set1_pd(70.000000),    
                                            _mm512_set1_pd(71.000000),    
                                            _mm512_set1_pd(72.000000),    
                                            _mm512_set1_pd(73.000000),    
                                            _mm512_set1_pd(74.000000),    
                                            _mm512_set1_pd(75.000000),    
                                            _mm512_set1_pd(76.000000),    
                                            _mm512_set1_pd(77.000000),    
                                            _mm512_set1_pd(78.000000),    
                                            _mm512_set1_pd(79.000000),    
                                            _mm512_set1_pd(80.000000),    
                                            _mm512_set1_pd(81.000000),    
                                            _mm512_set1_pd(82.000000),    
                                            _mm512_set1_pd(83.000000),    
                                            _mm512_set1_pd(84.000000),    
                                            _mm512_set1_pd(85.000000),    
                                            _mm512_set1_pd(86.000000),    
                                            _mm512_set1_pd(87.000000),    
                                            _mm512_set1_pd(88.000000),    
                                            _mm512_set1_pd(89.000000),    
                                            _mm512_set1_pd(90.000000),    
                                            _mm512_set1_pd(91.000000),    
                                            _mm512_set1_pd(92.000000),    
                                            _mm512_set1_pd(93.000000),    
                                            _mm512_set1_pd(94.000000),    
                                            _mm512_set1_pd(95.000000),    
                                            _mm512_set1_pd(96.000000),    
                                            _mm512_set1_pd(97.000000),    
                                            _mm512_set1_pd(98.000000),    
                                            _mm512_set1_pd(99.000000),    
                                            _mm512_set1_pd(100.00000)};   
	                const __m512d eps = _mm512_set1_pd(1.0e-12);
	                const __m512d pi  = _mm512_set1_pd(3.141592653589793e+00);
	                const __m512d C350= _mm512_set1_pd(3.5);
	                const __m512d C00 = _mm512_setzero_pd();
	                const __m512d C10 = _mm512_set1_pd(1.0);
	                const __m512d C20 = _mm512_set1_pd(2.0);
	                const __m512d C050= _mm512_set1_pd(0.5);
	                const __m512d C025= _mm512_set1_pd(0.25);
	                const __m512d C1772453850905516027298167483341 = 
	                                    _mm512_set1_pd(1.772453850905516027298167483341); // sqrt(PI)
	                register __m512d ei1,ei2,er,er0,er1,er2;
	                register __m512d eri,err,r,ss,w,w1,w2,x2,t0,t1;
	                register __m512d c0r,c0i,csr,csi;
	                int32_t k,n;
	                x2  = _mm512_mul_pd(zr,zr);
	                ctx = C10;
	                if(_mm512_cmp_pd_mask(zr,C350,_CMP_LE_OQ)) {
	                    er = C10;
	                    r  = er;
	                    for(k = 1; k != 100; ++k) {
	                        register __m512d ctx = LUT1[k];
	                        r   = _mm512_div_pd(_mm512_mul_pd(r,x2),ctx);
	                        er  = _mm512_add_pd(er,r);
	                        t0  = _mm512_abs_pd(_mm512_sub_pd(er,w));
	                        t1  = _mm512_mul_pd(eps,_mm512_abs_pd(er));
	                        if(_mm512_cmp_pd_mask(t0,t1,_CMP_LE_OQ)) {
	                           break;
	                        }
	                        w = er;
	                    }
	                    t0 = xexp(negate_zmm8r8(x2));
	                    t1 = _mm512_mul_pd(C1772453850905516027298167483341,x);
	                    c0r = _mm512_div_pd(C20,
	                                   _mm512_mul_pd(t0,t1));
	                    c0i = C00;
	                    er0= _mm512_mul_pd(c0r,er);                                 
	                }
	                else {
	                    er = C10;
	                    r  = er;
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),C05),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[0]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[1]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[2]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[3]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[4]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[5]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[6]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[7]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[8]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[9]),x2);
	                    er = _mm512_add_pd(er,r);
	                    r  = _mm512_div_pd(_mm512_mul_pd(
	                                             negate_zmm8r8(r),LUT1[10]),x2);
	                    er = _mm512_add_pd(er,r);
	                    t0 = xexp(negate_zmm8r8(x2));
	                    t1 = _mm512_mul_pd(x,
	                                   C1772453850905516027298167483341);
	                    c0r= _mm512_div_pd(t0,t1);
	                    c0i= C00;
	                    er0= _mm512_sub_pd(C10,_mm512_mul_pd(c0r,er));
	                }
	                
	                if(_mm512_cmp_pd_mask(y,C00,_CMP_EQ_OQ0)) {
	                    err = er0;
	                    eri = C00;
	                }
	                else {
	                    t0 = _mm512_mul_pd(C20,_mm512_mul_pd(x,y));
	                    csr= xcos(t0);
	                    ss = xsin(t0);
	                    t1 = _mm512_mul_pd(x,_mm512_mul_pd(C20,pi));
	                    t0 = xexp(negate_zmm8r8(x2));
	                    er1= _mm512_div_pd(_mm512_mul_pd(t0,
	                                                _mm512_sub_pd(C10,csr)),t1);
	                    ei1= _mm512_div_pd(_mm512_mul_pd(t0,ss),t1);
	                    er2= C00;
	                    for(n = 1; n != 100; ++n) {
	                        register __m512d xn = LUT2[n];
	                        t0                 = _mm512_mul_pd(xn,xn);
	                        t1                 = _mm512_mul_pd(xn,y);
	                        register __m512d zz1= xcosh(t1);
	                        register __m512d zz2= xsinh(t1);
	                        register __m512d zz3= xexp(negate_zmm8r8(
	                                                  _mm512_mul_pd(C025,t0)));
	                        register __m512d zz4= _mm512_fmadd_pd(xn,xn,
	                                                  _mm512_mul_pd(C40,x2));
	                        register __m512d zz5= _mm512_div_pd(zz3,zz4);
	                        register __m512d zz6= _mm512_fmsub_pd(C20,x,
	                                                       _mm512_mul_pd(C20,
	                                                                _mm512_mul_pd(x,zz1)));
	                        register __m512d zz7= _mm512_mul_pd(xn,
	                                                       _mm512_mul_pd(zz2,ss));
	                        register __m512d zz8= _mm512_fmadd_pd(zz5,zz6,zz7);
	                        er2                 = _mm512_add_pd(er2,zz8);
	                        register __m512d zz9= _mm512_abs_pd(
	                                             _mm512_div_pd(_mm512_sub_pd(er2,w1),er2));
	                        const __mmask8 m = _mm512_cmp_pd_mask(zz9,eps,_CMP_LT_OQ);
	                        if(m) {break;}
	                        w1 = er2;
	                    }
	                    c0r = _mm512_mul_pd(C20,_mm512_div_pd(xexp(negate_zmm8r8(x2))),pi);
	                    t0  = _mm512_add_pd(er0,er1);
	                    err = _mm512_fmadd_pd(er2,c0r,t0);
	                    ei2 = C00;
	                    for(n = 1; n != 100; ++n) {
	                        register __m512d xn = LUT2[n];
	                        t0                 = _mm512_mul_pd(xn,xn);
	                        t1                 = _mm512_mul_pd(xn,y);
	                        register __m512d zz1= xcosh(t1);
	                        register __m512d zz2= xsinh(t1);
	                        register __m512d zz3= xexp(negate_zmm8r8(
	                                                  _mm512_mul_pd(C025,t0)));
	                        register __m512d zz4= _mm512_fmadd_pd(xn,xn,
	                                                  _mm512_mul_pd(C40,x2));
	                        register __m512d zz5= _mm512_div_pd(zz3,zz4);
	                        register __m512d zz6= _mm512_mul_pd(_mm512_mul_pd(C20,x),
	                                                            _mm512_mul_pd(zz1,ss));
	                        register __m512d zz7= _mm512_mul_pd(xn,
	                                                        _mm512_mul_pd(zz2,csr));
	                        register __m512d zz8= _mm512_fmadd_pd(zz5,zz6,zz7);
	                        ei2                 = _mm512_add_pd(ei2,zz8);
	                        register __m512d zz9= _mm512_abs_pd(
	                                                 _mm512_div_pd(_mm512_sub_pd(ei2,w2),ei2));
	                        const __mmask8 m = _mm512_cmp_pd_mask(zz9,eps,_CMP_LT_OQ);
	                        if(m) {break;}
	                        w2 = ei2;
	                    }
	                    eri = _mm512_fmadd_pd(ei1,c0r,ei2);
	                }
	                t0    = _mm512_mul_pd(negate_zmm8r8(z),z);
	                t1    = xexp(t0);
	                *cerr = err;
	                *ceri = eri;
	                *cderr= _mm512_mul_pd(_mm512_set1_pd(1.128379167095512573896158903122),t1);
	                *cderi= C00;
	         }
	         
        
       } // math




} // gms































#endif /*__GMS_SPEC_FUNC_ZMM8R8_HPP__*/
