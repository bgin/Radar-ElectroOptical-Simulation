#ifndef __GMS_MATH_CONSTANTS_H__
#define __GMS_MATH_CONSTANTS_H__

 
#include <complex>

//#include <boost\math\special_functions\gamma.hpp>
namespace mathlib
{
	class MathConstants
	{
	public:

		//Value of PI - double precision.
		inline const static  double PI_DBL()
		{  
			//boost::math::tgamma<double>(0.5);
			return  3.1415926535897932e+00; //boost::math::constants::pi<double>();
		}

		//Value of PI - single precision.
		inline static  float  PI_FLT()
		{
			
			return  3.141592653e+00;//boost::math::constants::pi<float>();
		}

		inline static  double HALF_PI_DBL()
		{
			return (0.5 * PI_DBL());//boost::math::constants::half_pi<double>();
		}

		inline static  float HALF_PI_FLT()
		{
			return  0.5f * PI_FLT();//boost::math::constants::half_pi<float>();
		}

		inline static  double NEG_HALF_PI_DBL()
		{
			
			return -(HALF_PI_DBL());//(boost::math::constants::half_pi<double>());
		}

		inline static  float  NEG_HALF_PI_FLT()
		{
			return -(HALF_PI_FLT());//(boost::math::constants::half_pi<float>());
		}

		inline static  double NEG_PI_DBL()
		{
			return -(PI_DBL());//(boost::math::constants::pi<double>());
		}

		inline static  float NEG_PI_FLT()
		{
			return  -(PI_FLT());//(boost::math::constants::pi<float>());
		}

		inline static  double QUATR_PI_DBL()
		{
			return (0.25 * PI_DBL());
		}

		inline static  float QUATR_PI_FLT()
		{
			return (0.25f * PI_FLT());
		}

		inline static  double NEG_QUATR_PI_DBL()
		{
			return -(QUATR_PI_DBL());
		}

		inline static  float NEG_QUATR_PI_FLT()
		{
			return -(QUATR_PI_FLT());
		}

		inline static  double TWO_PI_DBL()
		{
			
			return  2.0 * PI_DBL();; //boost::math::constants::two_pi<double>();
		}

		inline static  float TWO_PI_FLT()
		{
			return 2.f * PI_FLT();//boost::math::constants::two_pi<float>();
		}

		inline static  double FOUR_PI_DBL()
		{
			return (4.0 * PI_DBL());
		}

		inline static  float  FOUR_PI_FLT()
		{
			return (4.0f * PI_FLT());
		}

		inline static  double INV_PI_DBL()
		{
			return (1.0 / PI_DBL());
		}

		inline static  float  INV_PI_FLT()
		{
			return (1.0f / PI_FLT());
		}

		inline static  double  INV_TWO_PI_DBL()
		{
			return (1.0 / TWO_PI_DBL());
		}

		inline static  float   INV_TWO_PI_FLT()
		{
			return (1.0f / INV_TWO_PI_FLT());
		}

		inline static  double  SQRT_PI_DBL()
		{
			
			return 1.7724538509055160e+00;//boost::math::constants::root_pi<double>();
		}

		inline static  float   SQRT_PI_FLT()
		{
			return 1.772453850e+00;//boost::math::constants::root_pi<float>();
		}

		inline static  double  INV_SQRT_PI_DBL()
		{
			return (1.0 / SQRT_PI_DBL());
		}

		inline static  float   INV_SQRT_PI_FLT()
		{
			return  (1.0f / SQRT_PI_FLT());
		}

		inline static std::complex<double>  SQRT_i_DBL()
		{
			std::complex<double> i = -1;
			return std::sqrt(i);
		}

		inline static  std::complex<float> SQRT_i_FLT()
		{
			std::complex<float> i = -1;
			return std::sqrt(i);
		}

		/*__forceinline static  double  SQRT_TWO_DBL()
		{
			
			return boost::math::constants::root_two<double>();
		}

		__forceinline static  float  SQRT_TWO_FLT()
		{
			return boost::math::constants::root_two<float>();
		}

		__forceinline static  double e_dbl()
		{
			return boost::math::constants::e<double>();
		}

		__forceinline static  float  e_flt()
		{
			return boost::math::constants::e<float>();
		}

		__forceinline static  double LN_TWO_DBL()
		{
			return boost::math::constants::ln_two<double>();
		}

		__forceinline static  float  LN_TWO_FLT()
		{
			return boost::math::constants::ln_two<float>();
		}*/
	};
}

#endif /*__GMS_MATH_CONSTANTS_H__*/
