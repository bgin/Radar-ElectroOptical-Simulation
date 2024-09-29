
#ifndef __GMS_SPEC_FUNCS_ZMM8R8_H__
#define __GMS_SPEC_FUNCS_ZMM8R8_H__ 180520230838


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

    const unsigned int GMS_SPEC_FUNCS_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_SPEC_FUNCS_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_SPEC_FUNCS_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_SPEC_FUNCS_ZMM8R8_FULLVER =
      1000U*GMS_SPEC_FUNCS_ZMM8R8_MAJOR+
      100U*GMS_SPEC_FUNCS_ZMM8R8_MINOR+
      10U*GMS_SPEC_FUNCS_ZMM8R8_MICRO;
    const char * const GMS_SPEC_FUNCS_ZMM8R8_CREATION_DATE = "18-05-2023 08:38 AM +00200 (THR 18 05 2023 GMT+2)";
    const char * const GMS_SPEC_FUNCS_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_SPEC_FUNCS_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_SPEC_FUNCS_ZMM8R8_DESCRIPTION   = "Special Functions library manually vectorized (AVX512 double)."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"

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
        
        
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besei0_zmm8r8(const __m512d);
	           
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besei0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	           
	           
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besei0_zmm8r8_u(const double * __restrict);
	           
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


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besei1_zmm8r8(const __m512d);
	           
	           
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besei1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	           
	           
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besei1_zmm8r8_u(const double * __restrict);
	           
	           
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besek0_zmm8r8(const __m512d);
	           
	           
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besek0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	           
	           
	            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besek0_zmm8r8_u(const double * __restrict);
	           
	           
	           
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


                   
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besek1_zmm8r8(const __m512d );
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besek1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	          
	          
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besek1_zmm8r8_u(const double * __restrict);
	           
	           
	           
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


                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besi0_zmm8r8(const __m512d);
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besi0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)); 
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besi0_zmm8r8_u(const double * __restrict);
	           
	           
	           
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besi1_zmm8r8(const __m512d);
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besi1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	         
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d besi1_zmm8r8_u(const double * __restrict); 
	           
	           
	           
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


                
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besj0_zmm8r8(const __m512d);
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besj0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	          
	          
	            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besj0_zmm8r8_u(const double * __restrict);
	           
	           
	           
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


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besj1_zmm8r8(const __m512d);
	          
	          
	             
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besj1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	          
	          
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besj1_zmm8r8_u(const double * __restrict); 
	           
	           
	           
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


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besk0_zmm8r8(const __m512d);
	          
	          
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besk0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)); 
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besk0_zmm8r8_u(const double * __restrict);
	           
	           
	           
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


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besk1_zmm8r8(const __m512d);
	          
	          
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besk1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	          
	          
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besk1_zmm8r8_u(const double * __restrict);
	           
	           
	           
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


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besy0_zmm8r8(const __m512d);
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besy0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besy0_zmm8r8_u(const double * __restrict);
	          
	          	    
	          	 
	          	           
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	       	   __m512d besy1_zmm8r8(const __m512d);
	          
	          
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besy1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64));
	          
	          
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d besy1_zmm8r8_u(const double * __restrict);
	          	          
	          	        
	          
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calcei_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calcei_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                   const int32_t);
	                                   
	                                   
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


                 
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calci0_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calci0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                 const int32_t);
	                                 
	                                 
	                                 
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calci1_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calci1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                   const int32_t);
	                                   
	                                   
	                                   
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
	               
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calck0_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	                                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calck0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                   const int32_t);
	                                   
	                                   
	                                              	                                                                                                                       	          
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


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calck1_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calck1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                 const int32_t);
	                                 
	                                 
	                                 
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calerf_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d calerf_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                   const int32_t);
	                                   
	                                   
	                                   
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


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d caljy0_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d caljy0_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                   const int32_t);
	                                   
	                                   
	                                   
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
	         
	         
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d caljy1_zmm8r8(const __m512d,
	                                 const int32_t);
	                                 
	                                 
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d caljy1_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                   const int32_t);  
	                                   
	                                   
	                                   
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
	          
	            
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          __m512d dawson_zmm8r8(const __m512d);	 
	          
	          
	          
	         __ATTR_HOT__
	         __ATTR_ALIGN__(32)
                 __ATTR_VECTORCALL__
	        __m512d dawson_zmm8r8_a(const double * 
	                   __restrict __ATTR_ALIGN__(64));
	                   
	                   
	                   
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
	        
	        
	         
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)
                __ATTR_VECTORCALL__
	        void cerf_zmm8r8(const __m512d,
	                         const __m512d,
	        __m512d * __restrict ,
	        __m512d * __restrict ,
	        __m512d * __restrict ,
	        __m512d * __restrict ); 
	        
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cerf_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ,
	                              const double * __restrict __ATTR_ALIGN__(64) ,
	                              double * __restrict __ATTR_ALIGN__(64) ,
	                              double * __restrict __ATTR_ALIGN__(64) ,
	                              double * __restrict __ATTR_ALIGN__(64) ,
	                              double * __restrict __ATTR_ALIGN__(64) );
	                              
	                              
	                              
	                              
/*
     !*****************************************************************************80
!
!! CERROR computes the error function for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
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
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CER, the function value.
*/	         
	         
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	            void cerror_zmm8r8(const __m512d,
	                              const __m512d,
	                              __m512d * __restrict,
	                              __m512d * __restrict);
	                              
	                              
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           void cerror_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64),
	                                const double * __restrict __ATTR_ALIGN__(64),
	                                double * __restrict __ATTR_ALIGN__(64),
	                                double * __restrict __ATTR_ALIGN__(64));
	                                
	                                
#if 0
  !*****************************************************************************80
!
!! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes the modified Bessel functions I0(z), I1(z), 
!    K0(z), K1(z), and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
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
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1, 
!    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z), 
!    and K1'(z).
#endif     


                 


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           void cik01_zmm8r8(const __m512d,
	                             const __m512d,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict,
	                             __m512 * __restrict);
	                             
	                             
	                                                                      
/*
 !*****************************************************************************80
!
!! CJY01: complex Bessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
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
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, 
!    CDY1, the values of J0(z), J0'(z), J1(z), J1'(z), Y0(z), Y0'(z), 
!    Y1(z), Y1'(z).
!
*/	      
	   
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	            void cjy01_zmm8r8(const __m512d,
	                             const __m512d,
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64),
	                             __m512d * __restrict __ATTR_ALIGN__(64));
	                             
	                             
	                                   
	                                                                                                                                                                                                        	           
/*
         
!*****************************************************************************80
!
!! CV0 computes the initial characteristic value of Mathieu functions.
!
!  Discussion:
!
!    This procedure computes the initial characteristic value of Mathieu 
!    functions for m <= 12 or q <= 300 or q <= m*m.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!   03 August 2012
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
!    Input, integer ( kind = 4 ) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer ( kind = 4 ) M, the order of the functions.
!
!    Input, real ( kind = 8 ) Q, the parameter of the functions.
!
!    Output, real ( kind = 8 ) A0, the characteristic value.
*/	


                
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          __m512d  cv0_case0_zmm8r8(const int32_t,
	                                   const int32_t,
	                                   const __m512d);
	                                   
	                                   
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cv0_case1_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);
	                                     
	                                     
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cv0_case2_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);
	                                     
	                                     
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cv0_case3_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);
	                                     
	                                     
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         __m512d  cv0_case4_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);      
	                                     
	                                     
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cv0_case5_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);
	                                     
	                                     
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          __m512d  cv0_case6_zmm8r8( const int32_t,
	                                     const int32_t,
	                                     const __m512d);   
	                                     
	                                     
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cv0_case7_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);
	                                     
	                                     
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cv0_case8_zmm8r8(const int32_t,
	                                     const int32_t,
	                                     const __m512d);   
	                                     
	                                     
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d  cv0_zmm8r8(const int32_t,
	                               const int32_t,
	                               const __m512d);     
	                               
	                               
/*
!*****************************************************************************80
!
!! CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
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
!    Input, integer ( kind = 4 ) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer ( kind = 4 ) M, the order of the Mathieu functions.
!
!    Input, real ( kind = 8 ) Q, the parameter value.
!
!    Output, real ( kind = 8 ) A0, the initial characteristic value.
*/


                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d cvql_zmm8r8(const int32_t,
	                               const int32_t,
	                               const __m512d);
	                               
	                               
	            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d cvql_zmm8r8_a(const int32_t,
	                                 const int32_t,
	                                 const double * __restrict __ATTR_ALIGN__(64));
	                                 
	                                 
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d cvql_zmm8r8_u(const int32_t,
	                                 const int32_t,
	                                 const double * __restrict);
	                                 
	                                 
/*
   !*****************************************************************************80
!
!! CVQM computes the characteristic value of Mathieu functions for q <= m*m.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
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
!    Input, integer ( kind = 4 ) M, the order of the Mathieu functions.
!
!    Input, real ( kind = 8 ) Q, the parameter value.
!
!    Output, real ( kind = 8 ) A0, the initial characteristic value.   
*/	


                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cvqm_zmm8r8(const int32_t,
	                                const __m512d);
	                                
	                                
	            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cvqm_zmm8r8_a(const int32_t,
	                                 const double * __restrict __ATTR_ALIGN__(64));	 
	                                 
	                                 
	            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           __m512d  cvqm_zmm8r8_u(const int32_t,
	                                 const double * __restrict);      
	                                 
	                                 
	                                                                                                                                                               
	                                   	          	           	           	           
	            
          
       } // math


} // gms


















#endif /*__GMS_SPEC_FUNCS_ZMM8R8_H__*/
