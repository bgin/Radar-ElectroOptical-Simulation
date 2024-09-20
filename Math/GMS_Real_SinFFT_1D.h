#ifndef __GMS_REAL_SINE_FFT_1D_H__
#define __GMS_REAL_SINE_FFT_1D_H__

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Real Sine FFT forward/backward Direction 1D class  - wrapper definition.
Wrapper around FFTPACK 5.1 library
@file Real_SinFFT_1D.h
@author: Bernard Gingold
@version:  1.0  26/10/2015
@description: in F77_FFT_DECLARATIONS.h
*/

#include <vector>
#include "GMS_fftpack51_declarations.h"

namespace   gms
{

	class RealSinFFT1D
	{

	public:


		/*
		@brief wrapper for F77 SINT1I routine.
		*/
		void      real_sinfft_init1D(int,
			std::vector<float> &,
			int,
			int);

		/*
		@brief wrapper for F77 SINT1B routine
		*/
		void       real_sinfft_backward1D(int,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			int);

		/*
		@brief wrapper for F77 SINT1F routine.
		*/
		void       real_sinfft_forward1D(int,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			std::vector<float> &,
			int,
			int);
			
	};

#include "GMS_Real_SinFFT_1D.inl"
}
#endif /*__GMS_REAL_SINE_FFT_1D_H__*/
