
#ifndef __GMS_FAST_MATH_FUNCS_H__
#define __GMS_FAST_MATH_FUNCS_H__

#include <math.h>
#include <type_traits>

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

    const unsigned int GMS_FAST_MATH_FUNCS_MAJOR = 1U;
    const unsigned int GMS_FAST_MATH_FUNCS_MINOR = 0U;
    const unsigned int GMS_FAST_MATH_FUNCS_MICRO = 0U;
    const unsigned int GMS_FAST_MATH_FUNCS_FULLVER =
      1000U*GMS_FAST_MATH_FUNCS_MAJOR+
      100U*GMS_FAST_MATH_FUNCS_MINOR+
      10U*GMS_FAST_MATH_FUNCS_MICRO;
    const char * const GMS_FAST_MATH_FUNCS_CREATION_DATE = "07-11-2023 14:40  +00200 (TUE 07 NOV 2023 GMT+2)";
    const char * const GMS_FAST_MATH_FUNCS_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_FAST_MATH_FUNCS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_FAST_MATH_FUNCS_DESCRIPTION   = "Simple Horner-scheme based approximations."

}

namespace gms {


   namespace math {
   
   
	template<typename T> class fast_math_funcs
	{
	public:

		/*template<bool B, typename T=void>
		using Enable_if = typename std::enable_if<B, T>::type;*/
			// Fast cotangent polynomial approximation.
		inline static  std::enable_if<std::is_floating_point<T>::value, T> &fastcot(const T);
		
		// Fast sine function polynomial approximation
		inline static  std::enable_if<std::is_floating_point<T>::value, T>  &fastsin(const T);

		// Fast cosine function polynomial approximation.
		inline static std::enable_if<std::is_floating_point<T>::value, T>  &fastcos(const T);
			

		// Fast tangent polynomial aaproximation
		inline 	static   std::enable_if<std::is_floating_point<T>::value, T>  &fasttan(const T);

		// Fast exponential function plynomial aaproximation.
		inline static std::enable_if<std::is_floating_point<T>::value, T>   &fastexp(const T);

		// Fast cosecant function polynomial aaproximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastcsc(const T);

		// Fast secant function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastsec(const T);

		// Fast arcsin function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastarcsin(const T);

		// Fast arctangent function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastarctan(const T);

		// Fast hyperbolic sine function polynomial aaproximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastsinh(const T);

		// Fast hyperbolic cosine function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastcosh(const T);

		// Fast hyperbolic tangent function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fasttanh(const T);

		// Fast hyperbolic cosecant function polynomial approximation
		inline static std::enable_if < std::is_floating_point<T>::value, T> &fastcsch(const T);

		// Fast hyperbolic secant function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastsech(const T);

		// Fast hyperblic contangent function polynomial approximation
		inline static std::enable_if<std::is_floating_point<T>::value, T> &fastcoth(const T);
		//::sin<
	

	    /* inline 	static constexpr bool isFloatingPoint()
		{
			return std::is_floating_point<T>::value;
		}*/
		/*static constexpr T TEST = 3.14;*/ //Weird bug and should be allowed
		
		//static constexpr auto TEST_PI = 3.14;
		
		
	
	};

#include "GMS_fast_math_funcs.inl"
    }
}
#endif /*__GMS_FAST_MATH_FUNCS_H__*/
