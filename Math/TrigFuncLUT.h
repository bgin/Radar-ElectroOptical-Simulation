#ifndef _TRIG_FUNC_LUT_H_
#define _TRIG_FUNC_LUT_H_

#include "MathLibDefinitions.h"


namespace    mathlib
{

	template<typename _Ty>  class SinCoeffLUT
	{

	public:

		
		
		/*
		@brief  static const array of 8*8 Sin(2*Pi*k*i/T) elements where T=8, k=1....8, i=0...7.
		*/
		static const   _Ty m_sSinCoeff8x8[64];

	
		/*
		@brief  static const array of 8*8 Sin(Pi*k*i/T) of half period elements where, T=8, k=1....8, i=0....7.
		*/
		static const  _Ty m_sSinCoeffH8x8[64];
		
		/*
		@brief  static const array of 16*16 Sin(2*Pi*k*i/T) elements where, T=16, k=1....16, i=0....15.
		*/
		static const  _Ty m_sSinCoeff16x16[256];

		/*
		@brief   static const array of 16*16 Sin(Pi*k*i/T) of half period elements  where, T=16, k=1....16, i=0....15.
		*/
		static const  _Ty m_sSinCoeffH16x16[256];

		/*
		@brief   static const array of 32*32 Sin(2*Pi*k*i/T) elements where, T=16, k=1....32, i=0....31.
		*/
		static const  _Ty m_sSinCoeff32x32[1024];

	};


	template<typename _Ty>  class CosCoeffLUT
	{

	public:


		/*
		@brief  static constexpr array of 8*8 Cos(2*Pi*k*i/T) where T=8, k=1....8, i=0...7.
		*/
		static const  _Ty m_sCosCoeff8x8[64];

		/*
		@brief  static const    array  of 8*8 Cos(Pi*k*i/T) half period elements where t=8, k=1....8, i=0....7.
		*/
		static const  _Ty m_sCosCoeffH8x8[64];

		/*
		@brief  static const    array  of 16*16 Cos(2*Pi*k*i/T) elements where, T=16, k=1....16, i=0....15.
		*/
		static const  _Ty m_sCosCoeff16x16[256];

		/*
		@brief   static const array of 16*16 Cos(Pi*k*i/T) of half period elements  where, T=16, k=1....16, i=0....15.
		*/
		static const  _Ty m_sCosCoeffH16x16[256];

		/*
		@brief   static const array of 32*32 Sin(2*Pi*k*i/T) elements where, T=16, k=1....32, i=0....31.
		*/
		static const  _Ty m_sCosCoeff32x32[1024];
	};

#include "TrigFuncLUT.inl"
}
#endif /*_TRIG_FUNC_LUT_H_*/